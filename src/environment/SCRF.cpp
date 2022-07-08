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

#include "SCRF.h"

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "chemistry/PhysicalConstants.h"
#include "chemistry/chemistry_utils.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/two_electron/ReactionPotential.h"
#include "scf_solver/KAIN.h"
#include "utils/print_utils.h"

#include <nlohmann/json.hpp>

using mrcpp::Printer;
using mrcpp::Timer;

using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;
using DerivativeOperator_p = std::shared_ptr<mrcpp::DerivativeOperator<3>>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {

SCRF::SCRF(Permittivity e,
           const Nuclei &N,
           PoissonOperator_p P,
           DerivativeOperator_p D,
           double orb_prec,
           int kain_hist,
           int max_iter,
           bool accelerate_Vr,
           std::string convergence_criterion,
           std::string algorithm,
           std::string density_type)
        : accelerate_Vr(accelerate_Vr)
        , convergence_criterion(convergence_criterion)
        , algorithm(algorithm)
        , density_type(density_type)
        , max_iter(max_iter)
        , history(kain_hist)
        , apply_prec(orb_prec)
        , conv_thrs(1.0)
        , mo_residual(-1.0)
        , epsilon(e)
        , rho_nuc(false)
        , rho_ext(false)
        , rho_tot(false)
        , Vr_n(false)
        , dVr_n(false)
        , Vr_nm1(false)
        , gamma_n(false)
        , dgamma_n(false)
        , gamma_nm1(false)
        , derivative(D)
        , poisson(P) {
    setDCavity();
    rho_nuc = chemistry::compute_nuclear_density(this->apply_prec, N, 1000);
}

SCRF::~SCRF() {
    mrcpp::clear(this->d_cavity, true);
}

void SCRF::clear() {
    this->rho_tot.free(NUMBER::Real);
}

double SCRF::setConvergenceThreshold(double prec) {
    // converge_thrs should be in the interval [prec, 1.0]
    this->conv_thrs = prec;
    bool dynamic = (this->convergence_criterion == "dynamic");
    if (dynamic and this->mo_residual > 10 * prec) this->conv_thrs = std::min(1.0, this->mo_residual);
    return this->conv_thrs;
}

void SCRF::setDCavity() {
    mrcpp::FunctionTree<3> *dx_cavity = new mrcpp::FunctionTree<3>(*MRA);
    mrcpp::FunctionTree<3> *dy_cavity = new mrcpp::FunctionTree<3>(*MRA);
    mrcpp::FunctionTree<3> *dz_cavity = new mrcpp::FunctionTree<3>(*MRA);
    d_cavity.push_back(std::make_tuple(1.0, dx_cavity));
    d_cavity.push_back(std::make_tuple(1.0, dy_cavity));
    d_cavity.push_back(std::make_tuple(1.0, dz_cavity));
    mrcpp::project<3>(this->apply_prec / 100, this->d_cavity, this->epsilon.getGradVector());
}

void SCRF::computeDensities(OrbitalVector &Phi) {
    Timer timer;
    resetQMFunction(this->rho_tot);
    Density rho_el(false);
    density::compute(this->apply_prec, rho_el, Phi, DensityType::Total);
    rho_el.rescale(-1.0);
    if (this->density_type == "electronic") {
        qmfunction::deep_copy(this->rho_tot, rho_el);
    } else if (this->density_type == "nuclear") {
        qmfunction::deep_copy(this->rho_tot, this->rho_nuc);
    } else {
        qmfunction::add(this->rho_tot, 1.0, rho_el, 1.0, this->rho_nuc, -1.0); // probably change this into a vector
    }
    print_utils::qmfunction(3, "Vacuum density", this->rho_tot, timer);
}

void SCRF::computeGamma(QMFunction &potential, QMFunction &out_gamma) {
    auto d_V = mrcpp::gradient(*derivative, potential.real());
    resetQMFunction(out_gamma);
    mrcpp::dot(this->apply_prec, out_gamma.real(), d_V, this->d_cavity);
    out_gamma.rescale(std::log((epsilon.getEpsIn() / epsilon.getEpsOut())) * (1.0 / (4.0 * mrcpp::pi)));
    mrcpp::clear(d_V, true);
}

QMFunction SCRF::solvePoissonEquation(const QMFunction &in_gamma) {
    QMFunction Poisson_func;
    QMFunction rho_eff;
    QMFunction first_term;
    QMFunction Vr;
    QMFunction eps_inv;
    eps_inv.alloc(NUMBER::Real);
    Vr.alloc(NUMBER::Real);

    this->epsilon.flipFunction(true);
    qmfunction::project(eps_inv, this->epsilon, NUMBER::Real, this->apply_prec / 100);
    this->epsilon.flipFunction(false);
    qmfunction::multiply(first_term, this->rho_tot, eps_inv, this->apply_prec);
    qmfunction::add(rho_eff, 1.0, first_term, -1.0, this->rho_tot, -1.0);

    qmfunction::add(Poisson_func, 1.0, in_gamma, 1.0, rho_eff, -1.0);
    mrcpp::apply(this->apply_prec, Vr.real(), *poisson, Poisson_func.real());

    return Vr;
}

void SCRF::accelerateConvergence(QMFunction &dfunc, QMFunction &func, KAIN &kain) {
    OrbitalVector phi_n(0);
    OrbitalVector dPhi_n(0);
    phi_n.push_back(Orbital(SPIN::Paired));
    dPhi_n.push_back(Orbital(SPIN::Paired));

    phi_n[0] = func;
    dPhi_n[0] = dfunc;

    kain.accelerate(this->apply_prec, phi_n, dPhi_n);

    func = phi_n[0];
    dfunc = dPhi_n[0]; // see if you can remove these two lines

    phi_n.clear();
    dPhi_n.clear();
}

void SCRF::nestedSCRF(QMFunction V_vac) {
    KAIN kain(this->history);
    kain.setLocalPrintLevel(10);

    mrcpp::print::separator(3, '-');

    double update = 10.0, norm = 1.0;
    int iter = 1;
    while (update >= this->conv_thrs && iter <= max_iter) {
        Timer t_iter;
        // solve the poisson equation
        QMFunction Vr_np1 = solvePoissonEquation(this->gamma_n);
        norm = Vr_np1.norm();

        // use a convergence accelerator
        resetQMFunction(this->dVr_n);
        qmfunction::add(this->dVr_n, 1.0, Vr_np1, -1.0, this->Vr_n, -1.0);
        update = dVr_n.norm();

        if (iter > 1 and this->history > 0 and this->accelerate_Vr) {
            accelerateConvergence(dVr_n, Vr_n, kain);
            Vr_np1.free(NUMBER::Real);
            qmfunction::add(Vr_np1, 1.0, Vr_n, 1.0, dVr_n, -1.0);
        }

        // set up for next iteration
        QMFunction V_tot;
        qmfunction::add(V_tot, 1.0, Vr_np1, 1.0, V_vac, -1.0);
        updateCurrentReactionPotential(Vr_np1); // push_back() maybe

        QMFunction gamma_np1;
        computeGamma(V_tot, gamma_np1);

        resetQMFunction(dgamma_n);
        qmfunction::add(this->dgamma_n, 1.0, gamma_np1, -1.0, this->gamma_n, -1.0);

        if (iter > 1 and this->history > 0 and (not this->accelerate_Vr)) {
            accelerateConvergence(dgamma_n, gamma_n, kain);
            gamma_np1.free(NUMBER::Real);
            qmfunction::add(gamma_np1, 1.0, gamma_n, 1.0, dgamma_n, -1.0);
        }

        updateCurrentGamma(gamma_np1);
        printConvergenceRow(iter, norm, update, t_iter.elapsed());
        iter++;
    }
    mrcpp::print::separator(3, '-');
    this->dVr_n.real().clear();
    this->dVr_n.real().setZero();
    this->dgamma_n.real().clear();
    this->dgamma_n.real().setZero();
    kain.clear();
}

void SCRF::printConvergenceRow(int i, double norm, double update, double time) const {
    auto pprec = Printer::getPrecision();
    auto w0 = Printer::getWidth() - 1;
    auto w1 = 9;
    auto w2 = 2 * w0 / 9;
    auto w3 = w0 - w1 - 2 * w2;

    std::string time_unit = " sec";
    if (time < 0.01) {
        time = time * 1000.0;
        time_unit = "  ms";
    } else if (time > 60.0) {
        time = time / 60.0;
        time_unit = " min";
    }

    std::stringstream o_txt;
    o_txt << " Iter " << std::setw(w1 - 6) << i;
    o_txt << " : " << std::setw(w3 - 3) << std::setprecision(2 * pprec) << std::scientific << norm;
    o_txt << std::setw(w2) << std::setprecision(pprec) << std::scientific << update;
    o_txt << std::setw(w2 - 4) << std::setprecision(2) << std::fixed << time << time_unit;
    println(3, o_txt.str());
}

QMFunction &SCRF::setup(double prec, const OrbitalVector_p &Phi) {
    this->apply_prec = prec;
    computeDensities(*Phi);
    Timer t_vac;
    QMFunction V_vac;
    V_vac.alloc(NUMBER::Real);
    mrcpp::apply(this->apply_prec, V_vac.real(), *poisson, this->rho_tot.real());
    print_utils::qmfunction(3, "Vacuum potential", V_vac, t_vac);

    // set up the zero-th iteration potential and gamma, so the first iteration gamma and potentials can be made

    Timer t_gamma;
    if ((not this->Vr_n.hasReal()) or (not this->gamma_n.hasReal())) {
        QMFunction gamma_0;
        QMFunction V_tot;
        computeGamma(V_vac, gamma_0);
        this->Vr_n = solvePoissonEquation(gamma_0);
        qmfunction::add(V_tot, 1.0, V_vac, 1.0, this->Vr_n, -1.0);
        computeGamma(V_tot, this->gamma_n);
    }

    // update the potential/gamma before doing anything with them

    if (accelerate_Vr) {
        QMFunction temp_Vr_n;
        qmfunction::add(temp_Vr_n, 1.0, this->Vr_n, 1.0, this->dVr_n, -1.0);
        qmfunction::deep_copy(this->Vr_n, temp_Vr_n);
        temp_Vr_n.free(NUMBER::Real);
        QMFunction V_tot;
        qmfunction::add(V_tot, 1.0, this->Vr_n, 1.0, V_vac, -1.0);
        resetQMFunction(this->gamma_n);
        computeGamma(V_tot, this->gamma_n);
    } else {
        QMFunction temp_gamma_n;
        qmfunction::add(temp_gamma_n, 1.0, this->gamma_n, 1.0, this->dgamma_n, -1.0);
        qmfunction::deep_copy(this->gamma_n, temp_gamma_n);
        temp_gamma_n.free(NUMBER::Real);
    }
    print_utils::qmfunction(3, "Initial gamma", this->gamma_n, t_gamma);

    Timer t_scrf;
    if (this->algorithm == "scrf") nestedSCRF(V_vac);
    print_utils::qmfunction(3, "Reaction potential", this->Vr_n, t_scrf);
    return this->Vr_n;
}

double SCRF::getNuclearEnergy() {
    return qmfunction::dot(this->rho_nuc, this->Vr_n).real();
}

double SCRF::getElectronicEnergy() {
    QMFunction rho_el;
    rho_el.alloc(NUMBER::Real);
    qmfunction::add(rho_el, 1.0, this->rho_tot, -1.0, this->rho_nuc, -1.0);
    return qmfunction::dot(rho_el, this->Vr_n).real();
}

double SCRF::getTotalEnergy() {
    return qmfunction::dot(this->rho_tot, this->Vr_n).real();
}

void SCRF::resetQMFunction(QMFunction &function) {
    if (function.hasReal()) function.free(NUMBER::Real);
    if (function.hasImag()) function.free(NUMBER::Imag);
    function.alloc(NUMBER::Real);
}

void SCRF::updateCurrentReactionPotential(QMFunction &Vr_np1) {
    resetQMFunction(this->Vr_nm1);
    qmfunction::deep_copy(this->Vr_nm1, this->Vr_n);
    resetQMFunction(this->Vr_n);
    qmfunction::deep_copy(this->Vr_n, Vr_np1);
    Vr_np1.free(NUMBER::Real);
}

void SCRF::updateCurrentGamma(QMFunction &gamma_np1) {
    resetQMFunction(this->gamma_nm1);
    qmfunction::deep_copy(this->gamma_nm1, this->gamma_n);
    resetQMFunction(this->gamma_n);
    qmfunction::deep_copy(this->gamma_n, gamma_np1);
    gamma_np1.free(NUMBER::Real);
}

void SCRF::printParameters() const {
    bool dynamic = (this->convergence_criterion == "dynamic");

    std::stringstream o_iter;
    if (this->max_iter > 0) {
        o_iter << this->max_iter;
    } else {
        o_iter << "Off";
    }

    std::stringstream o_kain;
    if (this->history > 0) {
        o_kain << this->history;
    } else {
        o_kain << "Off";
    }

    nlohmann::json data = {
        {"Method                ", "SCRF"},
        {"Max iterations        ", max_iter},
        {"KAIN solver           ", o_kain.str()},
        {"Density type          ", density_type},
        {"Dynamic threshold     ", (dynamic) ? "On" : "Off"},
    };

    mrcpp::print::separator(3, '~');
    print_utils::json(3, data, false);
    mrcpp::print::separator(3, '~');
}

} // namespace mrchem
