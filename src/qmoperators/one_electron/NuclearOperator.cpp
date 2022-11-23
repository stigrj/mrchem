/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "NuclearOperator.h"

#include<fstream>
#include<limits>
#include<nlohmann/json.hpp>

#include "analyticfunctions/NuclearFunction.h"
#include "chemistry/chemistry_utils.h"
#include "parallel.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/QMFunction.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/QMPotential.h"
#include "utils/print_utils.h"
#include "utils/math_utils.h"


using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
    /*/
    if (proj_charge >= 2) {
        // open file
        std::ifstream f("param_V.json");
        // skip lines starting from '#'
        while (f.peek() == '#') f.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
        // read json data from file
        nlohmann::json j = nlohmann::json::parse(f);
        // close file
        f.close();

        // print json contents
        //std::cout << j.dump() << std::endl;

        double sum = 0.0;
        for (auto item : j["element"]) {
            const auto name = static_cast<std::string>( item["name"] );
            const double epsilon = static_cast<double>( item["epsilon"] );
            const double rms = static_cast<double>( item["rms"] );
            std::cout << name << " | " << epsilon << " | " << rms << std::endl;
            sum += rms;
        }
        std::cout << "Sum rms: " << sum << std::endl;

        return 0;
    }
    */

NuclearOperator::NuclearOperator(const Nuclei &nucs, double proj_prec, double exponent, bool mpi_share, int proj_charge) {
    if (proj_charge == 1) {
        coulomb_HFYGB(nucs, proj_prec, mpi_share);
    } else if (proj_charge == 2) {
        projectFiniteNucleus(nucs, proj_prec, exponent, mpi_share);
    } else if (proj_charge == 3) {
        homogeneus_charge_sphere(nucs, proj_prec, mpi_share);
    } else if (proj_charge == 4) {
        gaussian(nucs, proj_prec, mpi_share);
    } else {
        applyFiniteNucleus(nucs, proj_prec, proj_prec, exponent, mpi_share);
    }
}

/** Compute finite nucleus potential by projecting analytic expression for the potential.
*   Works only for single nucleus for now!
*/
void NuclearOperator::coulomb_HFYGB(const Nuclei &nucs, double proj_prec, bool mpi_share) {
    Timer t_tot;
    mrcpp::print::header(0, "Projecting nuclear potential coulomb_HFYGB");
 
    auto f_smooth = [&nucs,proj_prec] (const mrcpp::Coord<3> &r) -> double {
        double tmp_sum = 0.0;
        //for (int k = 0; k < nucs.size(); k++)
        for (auto &nuc : nucs) {
            const auto Z_I = nuc.getCharge();
            const auto &R_I = nuc.getCoord();
            const double R = math_utils::calc_distance(r, R_I);
            //double factor = (0.00435 * proj_prec / Z_I**5)**(1./3.);
            const double Z_I_2 = Z_I * Z_I;
            const double Z_I_5 = Z_I_2 * Z_I_2 * Z_I;
            const double factor = std::cbrt(0.00435 * proj_prec / Z_I_5);
            //double u = std::erf(R)/R + (1/(3*std::sqrt(mrcpp::pi,0.5)))*(std::exp(-(R**2)) + 16*std::exp(-4*R**2));
            const double u = std::erf(R)/R + (1/(3*std::sqrt(mrcpp::pi)))*(std::exp(-(R*R)) + 16*std::exp(-4*R*R));
            tmp_sum += -(Z_I * u) / factor;
        }
        return tmp_sum;
    };


    Timer t_pot;
    auto V_nuc = std::make_shared<QMPotential>(1, mpi_share);
    qmfunction::project(*V_nuc, f_smooth, NUMBER::Real, proj_prec);
    t_pot.stop();

    // Invoke operator= to assign *this operator
    RankZeroOperator &O = (*this);
    O = V_nuc;
    O.name() = "V_nuc";
}



void NuclearOperator::projectFiniteNucleus(const Nuclei &nucs, double proj_prec, double exponent, bool mpi_share) {
    Timer t_tot;
    mrcpp::print::header(2, "Projecting nuclear potential");

    std::ifstream f("param_V.json");
    while (f.peek() == '#') f.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    nlohmann::json param_V = nlohmann::json::parse(f);
    f.close();

    auto f_smooth = [&nucs, &param_V] (const mrcpp::Coord<3> &r) -> double {
        double tmp_sum = 0.0;
        //for (int k = 0; k < nucs.size(); k++)
        //const auto A = nucs.getElement(); 
        for (auto &nuc : nucs) {
            const auto &A = nuc.getElement().getSymbol(); // or getName() ?
            // Here I want that if A is equal item element in json file than take the value expoent for such element.
            // exponent = exponentjson // for the element X = A
            double exponent = 0.0;
            for (auto item : param_V["element"]) {
                const auto name = static_cast<std::string>( item["name"] );
                if (name == A) {
                    exponent = static_cast<double>( item["xi"] ); // is it exponent?
                    break;
                }
            }
            const double sqrt_exp = std::sqrt(exponent);
            const auto Z_I = nuc.getCharge();
            const auto &R_I = nuc.getCoord();
            const double R = math_utils::calc_distance(r, R_I);
            tmp_sum += -(Z_I / R) * std::erf(sqrt_exp * R);
        }
        return tmp_sum;
    };

    Timer t_pot;
    auto V_nuc = std::make_shared<QMPotential>(1, mpi_share);
    qmfunction::project(*V_nuc, f_smooth, NUMBER::Real, proj_prec);
    t_pot.stop();

    // Invoke operator= to assign *this operator
    RankZeroOperator &O = (*this);
    O = V_nuc;
    O.name() = "V_nuc";
}



void NuclearOperator::homogeneus_charge_sphere(const Nuclei &nucs, double proj_prec, bool mpi_share) {
    Timer t_tot;
    mrcpp::print::header(0, "Projecting nuclear potential homogeneus_charge_sphere");

    std::ifstream f("param_V.json");
    while (f.peek() == '#') f.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    nlohmann::json param_V = nlohmann::json::parse(f);
    f.close();
 
    auto f_smooth = [&nucs,&param_V] (const mrcpp::Coord<3> &r) -> double {
        double tmp_sum = 0.0;
        //for (int k = 0; k < nucs.size(); k++) 
        //const auto A = nucs.getElement();
        //read json file
        for (auto &nuc : nucs) {
            const auto &A = nuc.getElement().getSymbol(); // or getName() ?
            //const RMS = RMSjson //for the elelemnt X = A; 
            double RMS = 0.0;
            for (auto item : param_V["element"]) {
                const auto name = static_cast<std::string>( item["name"] );
                if (name == A) {
                    RMS = static_cast<double>( item["rms"] );
                    break;
                }
            }
            const auto RMS2 = RMS*RMS;
            const auto Z_I = nuc.getCharge();
            const auto &R_I = nuc.getCoord();
            double R = math_utils::calc_distance(r, R_I);
            double R0 = std::sqrt(RMS2*(5.0/3.0));
            double prec, factor;
            if (R <= R0) {
                prec = -Z_I / (2.0*R0);
                factor = 3.0 - (R*R)/(R0*R0);
            } else { 
                prec = -Z_I / R;
                factor = 1.0;
            }
            tmp_sum += prec / factor;
        }
        return tmp_sum;
    };


    Timer t_pot;
    auto V_nuc = std::make_shared<QMPotential>(1, mpi_share);
    qmfunction::project(*V_nuc, f_smooth, NUMBER::Real, proj_prec);
    t_pot.stop();

    // Invoke operator= to assign *this operator
    RankZeroOperator &O = (*this);
    O = V_nuc;
    O.name() = "V_nuc";
}


void NuclearOperator::gaussian(const Nuclei &nucs, double proj_prec, bool mpi_share) {
    Timer t_tot;
    mrcpp::print::header(0, "Projecting nuclear potential homogeneus_charge_sphere");

    std::ifstream f("param_V.json");
    while (f.peek() == '#') f.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
    nlohmann::json param_V = nlohmann::json::parse(f);
    f.close();
 
    auto f_smooth = [&nucs,&param_V] (const mrcpp::Coord<3> &r) -> double {
        double tmp_sum = 0.0;
        //for (int k = 0; k < nucs.size(); k++) 
        for (auto &nuc : nucs) {
            //const auto A = nucs.getElement();
            const auto &A = nuc.getElement().getSymbol(); // or getName() ?
            //read json file
            //const epsilon = epsilonjson //for the elelemnt X = A;
            double epsilon = 0.0;
            for (auto item : param_V["element"]) {
                const auto name = static_cast<std::string>( item["name"] );
                if (name == A) {
                    epsilon = static_cast<double>( item["epsilon"] ); // is it exponent?
                    break;
                }
            }
            const auto Z_I = nuc.getCharge();
            const auto &R_I = nuc.getCoord();
            const double R = math_utils::calc_distance(r, R_I);
            const double prec = -Z_I / R;
            const double u = std::erf(std::sqrt(epsilon) * R);
            tmp_sum += prec * u;
        }
        return tmp_sum;
    };


    Timer t_pot;
    auto V_nuc = std::make_shared<QMPotential>(1, mpi_share);
    qmfunction::project(*V_nuc, f_smooth, NUMBER::Real, proj_prec);
    t_pot.stop();

    // Invoke operator= to assign *this operator
    RankZeroOperator &O = (*this);
    O = V_nuc;
    O.name() = "V_nuc";
}


/** Compute finite nucleus potential by projecting Gaussian charge distribution and then apply Poisson.
*   Should work for any number of nuclei, but all get the same exponent.
*/
void NuclearOperator::applyFiniteNucleus(const Nuclei &nucs, double proj_prec, double apply_prec, double exponent, bool mpi_share) {
    Timer t_tot;
    mrcpp::print::header(2, "Projecting nuclear density");
    println(0, "Nuclear exponent: " << exponent);

    mrcpp::PoissonOperator P(*MRA, proj_prec);

    Timer t_rho;
    auto rho_nuc = chemistry::compute_nuclear_density(proj_prec, nucs, exponent);
    rho_nuc.rescale(-1.0);
    t_rho.stop();

    Timer t_pot;
    auto V_nuc = std::make_shared<QMPotential>(1, mpi_share);
    V_nuc->alloc(NUMBER::Real);
    mrcpp::apply(apply_prec, V_nuc->real(), P, rho_nuc.real());
    t_pot.stop();

    t_tot.stop();
    print_utils::qmfunction(2, "Apply Poisson", *V_nuc, t_pot);
    mrcpp::print::footer(2, t_tot, 2);

    // Invoke operator= to assign *this operator
    RankZeroOperator &O = (*this);
    O = V_nuc;
    O.name() = "V_nuc";
}

/*! @brief NuclearOperator represents the function: sum_i Z_i/|r - R_i|
 *  @param nucs: Collection of nuclei that defines the potential
 *  @param proj_prec: Precision for projection of analytic function
 *  @param smooth_prec: Precision for smoothing of analytic function
 *  @param mpi_share: Should MPI ranks on the same machine share this function?
 */
NuclearOperator::NuclearOperator(const Nuclei &nucs, double proj_prec, double smooth_prec, bool mpi_share) {
    if (proj_prec < 0.0) MSG_ABORT("Negative projection precision");
    if (smooth_prec < 0.0) smooth_prec = proj_prec;

    Timer t_tot;
    mrcpp::print::header(2, "Projecting nuclear potential");

    // Setup local analytic function
    Timer t_loc;
    NuclearFunction f_loc;
    setupLocalPotential(f_loc, nucs, smooth_prec);

    // Scale precision by charge, since norm of potential is ~ to charge
    double Z_tot = 1.0 * chemistry::get_total_charge(nucs);
    double Z_loc = 1.0 * chemistry::get_total_charge(f_loc.getNuclei());
    double tot_prec = proj_prec / std::min(1.0 * Z_tot, std::sqrt(2.0 * Z_tot));
    double loc_prec = proj_prec / std::max(1.0, Z_loc); // relative prec

    // Scale precision by box size, so that accuracy is independent of box size
    double vol = 1.0;
    for (int i = 0; i < 3; i++) vol *= MRA->getWorldBox().getBoxLength(i);
    vol /= 262144;                   // we use as reference a cube 64x64x64
    vol = std::max(1.0, vol);        // do not scale for smaller boxes
    loc_prec /= pow(vol, 1.0 / 6.0); // norm of 1/r over the box ~ root_6(Volume)

    // Project local potential
    QMFunction V_loc(false);
    qmfunction::project(V_loc, f_loc, NUMBER::Real, loc_prec);
    t_loc.stop();
    mrcpp::print::separator(2, '-');
    print_utils::qmfunction(2, "Local potential", V_loc, t_loc);

    // Collect local potentials
    Timer t_com;
    auto V_tot = std::make_shared<QMPotential>(1, mpi_share);
    allreducePotential(tot_prec, *V_tot, V_loc);
    t_com.stop();

    t_tot.stop();
    print_utils::qmfunction(2, "Allreduce potential", *V_tot, t_com);
    mrcpp::print::footer(2, t_tot, 2);

    // Invoke operator= to assign *this operator
    RankZeroOperator &O = (*this);
    O = V_tot;
    O.name() = "V_nuc";
}

void NuclearOperator::setupLocalPotential(NuclearFunction &f_loc, const Nuclei &nucs, double smooth_prec) const {
    int pprec = Printer::getPrecision();
    int w0 = Printer::getWidth() - 1;
    int w1 = 5;
    int w2 = 8;
    int w3 = 2 * w0 / 9;
    int w4 = w0 - w1 - w2 - 3 * w3;

    std::stringstream o_head;
    o_head << std::setw(w1) << "N";
    o_head << std::setw(w2) << "Atom";
    o_head << std::string(w4, ' ');
    o_head << std::setw(w3) << "Charge";
    o_head << std::setw(w3) << "Precision";
    o_head << std::setw(w3) << "Smoothing";

    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    for (int k = 0; k < nucs.size(); k++) {
        const Nucleus &nuc = nucs[k];
        double Z = nuc.getCharge();
        double c = detail::nuclear_potential_smoothing(smooth_prec, Z);

        // All projection must be done on grand master in order to be exact
        int proj_rank = (mpi::numerically_exact) ? 0 : k % mpi::orb_size;
        if (mpi::orb_rank == proj_rank) f_loc.push_back(nuc, c);

        std::stringstream o_row;
        o_row << std::setw(w1) << k;
        o_row << std::setw(w2) << nuc.getElement().getSymbol();
        o_row << std::string(w4, ' ');
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << Z;
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << smooth_prec;
        o_row << std::setw(w3) << std::setprecision(pprec) << std::scientific << c;
        println(2, o_row.str());
    }
}

void NuclearOperator::allreducePotential(double prec, QMFunction &V_tot, QMFunction &V_loc) const {
    // Add up local contributions into the grand master
    mpi::reduce_function(prec, V_loc, mpi::comm_orb);
    if (mpi::grand_master()) {
        // If numerically exact the grid is huge at this point
        if (mpi::numerically_exact) V_loc.crop(prec);
    }

    if (not V_tot.hasReal()) V_tot.alloc(NUMBER::Real);
    if (V_tot.isShared()) {
        int tag = 3141;
        // MPI grand master distributes to shared masters
        mpi::broadcast_function(V_loc, mpi::comm_sh_group);
        if (mpi::share_master()) {
            // MPI shared masters copies the function into final memory
            mrcpp::copy_grid(V_tot.real(), V_loc.real());
            mrcpp::copy_func(V_tot.real(), V_loc.real());
        }
        // MPI share masters distributes to their sharing ranks
        mpi::share_function(V_tot, 0, tag, mpi::comm_share);
    } else {
        // MPI grand master distributes to all ranks
        mpi::broadcast_function(V_loc, mpi::comm_orb);
        // All MPI ranks copies the function into final memory
        mrcpp::copy_grid(V_tot.real(), V_loc.real());
        mrcpp::copy_func(V_tot.real(), V_loc.real());
    }
}

} // namespace mrchem
