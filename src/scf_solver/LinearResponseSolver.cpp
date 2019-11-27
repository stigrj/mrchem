/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "HelmholtzVector.h"
#include "KAIN.h"
#include "LinearResponseSolver.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

#include "qmoperators/two_electron/FockOperator.h"
#include "utils/print_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

LinearResponseSolver::LinearResponseSolver(bool dyn, FockOperator &F_0, OrbitalVector &Phi_0, ComplexMatrix &F_mat_0)
        : dynamic(dyn)
        , f_oper_0(&F_0)
        , f_mat_0(&F_mat_0)
        , phi_0(&Phi_0) {}

/** @brief Run orbital optimization
 *
 * Optimize orbitals until convergence thresholds are met. This algorithm iterates
 * the Sternheimer response equations in integral form. Common implementation for
 * static and dynamic response. Main points of the algorithm:
 *
 * Pre SCF: setup Helmholtz operators with unperturbed energies
 *
 *  1) Setup perturbed Fock operator
 *  2) For X and Y orbitals do:
 *     a) Apply Helmholtz operator on all orbitals
 *     b) Project out occupied space (1 - rho_0)
 *     c) Compute updates and errors
 *     d) Compute KAIN updates
 *  4) Compute property
 *  5) Check for convergence
 *
 */
bool LinearResponseSolver::optimize(double omega, FockOperator &F_1, OrbitalVector &X_n, OrbitalVector &Y_n) {
    printParameters(omega, F_1.perturbation().name());

    // Setup KAIN accelerators
    KAIN kain_x(this->history);
    KAIN kain_y(this->history);
    OrbitalVector &Phi_0 = *this->phi_0;
    ComplexMatrix &F_mat_0 = *this->f_mat_0;
    ComplexMatrix F_mat_x = F_mat_0 + omega * ComplexMatrix::Identity(Phi_0.size(), Phi_0.size());
    ComplexMatrix F_mat_y = F_mat_0 - omega * ComplexMatrix::Identity(Phi_0.size(), Phi_0.size());

    RankZeroTensorOperator V_0 = this->f_oper_0->potential();
    RankZeroTensorOperator V_1 = F_1.potential() + F_1.perturbation();

    double err_o = 1.0;
    double err_t = 1.0;
    DoubleVector errors_x = DoubleVector::Zero(Phi_0.size());
    DoubleVector errors_y = DoubleVector::Zero(Phi_0.size());

    this->error.push_back(err_t);
    this->property.push_back(0.0);

    // Setup Helmholtz operators (fixed, based on unperturbed system)
    double helm_prec = getHelmholtzPrec();
    HelmholtzVector H_x(helm_prec, F_mat_x.real().diagonal());
    HelmholtzVector H_y(helm_prec, F_mat_y.real().diagonal());
    ComplexMatrix L_mat_x = H_x.getLambdaMatrix();
    ComplexMatrix L_mat_y = H_y.getLambdaMatrix();

    auto plevel = Printer::getPrintLevel();
    if (plevel < 1) {
        printConvergenceHeader("Symmetric property");
        printConvergenceRow(0);
    }

    int nIter = 0;
    bool converged = false;
    while (nIter++ < this->maxIter or this->maxIter < 0) {
        std::stringstream o_header;
        o_header << "SCF cycle " << nIter;
        mrcpp::print::header(1, o_header.str(), 0, '#');
        mrcpp::print::separator(2, ' ', 1);

        // Initialize SCF cycle
        Timer t_scf, t_lap;
        double orb_prec = adjustPrecision(err_o);

        // Setup perturbed Fock operator (including V_1)
        F_1.setup(orb_prec);

        if (dynamic and plevel == 1) mrcpp::print::separator(1, '-');

        { // Iterate X orbitals
            // Compute argument: psi_i = sum_j [L-F]_ij*x_j + (1 - rho_0)V_1(phi_i)
            Timer t_arg;
            mrcpp::print::header(2, "Computing Helmholtz argument");
            t_lap.start();
            OrbitalVector Psi_1 = V_1(Phi_0);
            mrcpp::print::time(2, "Applying V_1", t_lap);

            t_lap.start();
            orbital::orthogonalize(this->orth_prec, Psi_1, Phi_0);
            mrcpp::print::time(2, "Projecting (1 - rho_0)", t_lap);

            t_lap.start();
            OrbitalVector Psi_2 = orbital::rotate(L_mat_x - F_mat_x, Phi_0);
            mrcpp::print::time(2, "Rotating orbitals", t_lap);

            OrbitalVector Psi = orbital::add(1.0, Psi_1, 1.0, Psi_2, -1.0);
            Psi_1.clear();
            Psi_2.clear();
            mrcpp::print::footer(2, t_arg, 2);
            if (plevel == 1) mrcpp::print::time(1, "Computing Helmholtz argument", t_arg);

            // Apply Helmholtz operators
            OrbitalVector X_np1 = H_x.apply(V_0, X_n, Psi);
            Psi.clear();

            // Projecting (1 - rho_0)X
            mrcpp::print::header(2, "Projecting occupied space");
            t_lap.start();
            orbital::orthogonalize(this->orth_prec, X_np1, Phi_0);
            mrcpp::print::time(2, "Projecting (1 - rho_0)", t_lap);
            mrcpp::print::footer(2, t_lap, 2);
            if (plevel == 1) mrcpp::print::time(1, "Projecting occupied space", t_lap);

            // Compute update and errors
            OrbitalVector dX_n = orbital::add(1.0, X_np1, -1.0, X_n);
            errors_x = orbital::get_norms(dX_n);
            X_np1.clear();

            // Compute KAIN update:
            kain_x.accelerate(orb_prec, X_n, dX_n);

            // Prepare for next iteration
            X_n = orbital::add(1.0, X_n, 1.0, dX_n);

            DoubleMatrix norms = DoubleMatrix::Zero(X_n.size(), 3);
            for (int n = 0; n < X_n.size(); n++) {
                norms(n, 0) = X_n[n].norm();
                if (X_n[n].hasReal()) norms(n, 1) = std::sqrt(X_n[n].real().getSquareNorm());
                if (X_n[n].hasImag()) norms(n, 2) = std::sqrt(X_n[n].imag().getSquareNorm());
            }
            mpi::allreduce_matrix(norms, mpi::comm_orb);
            println(0, norms);
        }

        if (dynamic and plevel == 1) mrcpp::print::separator(1, '-');

        if (dynamic) { // Iterate Y orbitals
            // Compute argument: psi_i = sum_j [L-F]_ij*y_j + (1 - rho_0)V_1.dagger(phi_i)
            Timer t_arg;
            mrcpp::print::header(2, "Computing Helmholtz argument");
            t_lap.start();
            OrbitalVector Psi_1 = V_1.dagger(Phi_0);
            mrcpp::print::time(2, "Applying V_1.dagger()", t_lap);

            t_lap.start();
            orbital::orthogonalize(this->orth_prec, Psi_1, Phi_0);
            mrcpp::print::time(2, "Projecting (1 - rho_0)", t_lap);

            t_lap.start();
            OrbitalVector Psi_2 = orbital::rotate(L_mat_y - F_mat_y, Phi_0);
            mrcpp::print::time(2, "Rotating orbitals", t_lap);

            OrbitalVector Psi = orbital::add(1.0, Psi_1, 1.0, Psi_2, -1.0);
            Psi_1.clear();
            Psi_2.clear();
            mrcpp::print::footer(2, t_arg, 2);
            if (plevel == 1) mrcpp::print::time(1, "Computing Helmholtz argument", t_arg);

            // Apply Helmholtz operators
            OrbitalVector Y_np1 = H_y.apply(V_0, Y_n, Psi);
            Psi.clear();

            // Projecting (1 - rho_0)X
            mrcpp::print::header(2, "Projecting occupied space");
            t_lap.start();
            orbital::orthogonalize(this->orth_prec, Y_np1, Phi_0);
            mrcpp::print::time(2, "Projecting (1 - rho_0)", t_lap);
            mrcpp::print::footer(2, t_lap, 2);
            if (plevel == 1) mrcpp::print::time(1, "Projecting occupied space", t_lap);

            // Compute update and errors
            OrbitalVector dY_n = orbital::add(1.0, Y_np1, -1.0, Y_n);
            errors_y = orbital::get_norms(dY_n);
            Y_np1.clear();

            // Compute KAIN update:
            kain_y.accelerate(orb_prec, Y_n, dY_n);

            // Prepare for next iteration
            Y_n = orbital::add(1.0, Y_n, 1.0, dY_n);

            DoubleMatrix norms = DoubleMatrix::Zero(Y_n.size(), 3);
            for (int n = 0; n < Y_n.size(); n++) {
                norms(n, 0) = Y_n[n].norm();
                if (Y_n[n].hasReal()) norms(n, 1) = std::sqrt(Y_n[n].real().getSquareNorm());
                if (Y_n[n].hasImag()) norms(n, 2) = std::sqrt(Y_n[n].imag().getSquareNorm());
            }
            mpi::allreduce_matrix(norms, mpi::comm_orb);
            println(0, norms);
        }

        // Compute property
        mrcpp::print::header(2, "Computing symmetric property");
        t_lap.start();
        double prop = F_1.perturbation().trace(Phi_0, X_n, Y_n).real();
        this->property.push_back(prop);
        mrcpp::print::footer(2, t_lap, 2);
        if (plevel == 1) mrcpp::print::time(1, "Computing symmetric property", t_lap);

        // Clear perturbed Fock operator
        F_1.clear();

        // Compute errors
        auto err_p = std::abs(getUpdate(this->property, nIter, false));
        err_o = std::max(errors_x.maxCoeff(), errors_y.maxCoeff());
        err_t = std::sqrt(errors_x.dot(errors_x) + errors_y.dot(errors_y));
        converged = checkConvergence(err_o, err_p);
        this->error.push_back(err_t);

        // Finalize SCF cycle
        if (plevel < 1) printConvergenceRow(nIter);
        printOrbitals(orbital::get_norms(X_n), errors_x, X_n, 1);
        if (dynamic) printOrbitals(orbital::get_norms(Y_n), errors_y, Y_n, 1, false);
        mrcpp::print::separator(1, '-');
        printResidual(err_t, converged);
        mrcpp::print::separator(2, '=', 2);
        printProperty();
        printMemory();

        mrcpp::print::footer(1, t_scf, 2, '#');
        mrcpp::print::separator(2, ' ', 2);

        if (converged) break;
    }

    printConvergence(converged, "Symmetric property");
    reset();

    return converged;
}

/** @brief Pretty printing of the computed property with update */
void LinearResponseSolver::printProperty() const {
    double prop_0(0.0), prop_1(0.0);
    int iter = this->property.size();
    if (iter > 1) prop_0 = this->property[iter - 2];
    if (iter > 0) prop_1 = this->property[iter - 1];

    int w0 = (Printer::getWidth() - 1);
    int w1 = 20;
    int w2 = w0 / 3;
    int w3 = 8;
    int w4 = w0 - w1 - w2 - w3;

    std::stringstream o_head;
    o_head << std::setw(w1) << " ";
    o_head << std::setw(w2) << "Value";
    o_head << std::setw(w4) << "Update";
    o_head << std::setw(w3) << "Done";

    mrcpp::print::separator(2, '=');
    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    printUpdate(1, " Symmetric property", prop_1, prop_1 - prop_0, this->propThrs);
    mrcpp::print::separator(2, '=', 2);
}

void LinearResponseSolver::printParameters(double omega, const std::string &oper) const {
    std::stringstream o_calc;
    if (this->dynamic) {
        o_calc << "Dynamic linear response";
    } else {
        o_calc << "Static linear response";
    }
    std::stringstream o_loc;
    if (this->localize) {
        o_loc << "On";
    } else {
        o_loc << "Off";
    }

    std::stringstream o_omega;
    o_omega << std::setprecision(5) << std::fixed << omega << " (au)";

    std::stringstream o_kain;
    if (this->history > 0) {
        o_kain << this->history;
    } else {
        o_kain << "Off";
    }
    std::stringstream o_iter;
    if (this->maxIter > 0) {
        o_iter << this->maxIter;
    } else {
        o_iter << "Off";
    }

    std::stringstream o_prec_0, o_prec_1;
    o_prec_0 << std::setprecision(5) << std::scientific << this->orbPrec[0];
    o_prec_1 << std::setprecision(5) << std::scientific << this->orbPrec[1];

    std::stringstream o_thrs_p;
    if (this->propThrs < 0.0) {
        o_thrs_p << "Off";
    } else {
        o_thrs_p << std::setprecision(5) << std::scientific << this->propThrs;
    }

    std::stringstream o_thrs_o;
    if (this->orbThrs < 0.0) {
        o_thrs_o << "Off";
    } else {
        o_thrs_o << std::setprecision(5) << std::scientific << this->orbThrs;
    }

    std::stringstream o_helm;
    if (this->helmPrec < 0.0) {
        o_helm << "Dynamic";
    } else {
        o_helm << std::setprecision(5) << std::scientific << this->helmPrec;
    }

    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation        ", o_calc.str());
    if (dynamic) print_utils::text(0, "Frequency          ", o_omega.str());
    print_utils::text(0, "Method             ", this->methodName);
    print_utils::text(0, "Perturbation       ", oper);
    print_utils::text(0, "Max iterations     ", o_iter.str());
    print_utils::text(0, "KAIN solver        ", o_kain.str());
    print_utils::text(0, "Localization       ", o_loc.str());
    print_utils::text(0, "Start precision    ", o_prec_0.str());
    print_utils::text(0, "Final precision    ", o_prec_1.str());
    print_utils::text(0, "Helmholtz precision", o_helm.str());
    print_utils::text(0, "Property threshold ", o_thrs_p.str());
    print_utils::text(0, "Orbital threshold  ", o_thrs_o.str());
    mrcpp::print::separator(0, '~', 2);
}

} // namespace mrchem
