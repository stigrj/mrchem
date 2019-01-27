/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2019 Stig Rune Jensen, Jonas Juselius, Luca Frediani, and contributors.
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

#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"

#include "HelmholtzVector.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/RankZeroTensorOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

/** @brief HelmholtzVector constructor
 *
 * This will set the build precision of the Helmholtz operators and the vector
 * of lambda parameters that will be used in the subsequent application. No
 * operators are constructed at this point, they are produced on-the-fly in
 * the application.
 */
HelmholtzVector::HelmholtzVector(double pr, const DoubleVector &l)
        : prec(pr) {
    this->lambda = l;
    for (int i = 0; i < this->lambda.size(); i++) {
        if (this->lambda(i) > 0.0) this->lambda(i) = -0.5;
    }
}

/** @brief Apply Helmholtz operator component wise on OrbitalVector
 *
 * This will construct a separate Helmholtz operator for each of the entries
 * in the OrbitalVector based on the corresponding lambda_i parameter in the
 * HelmholtzVector. Computes output as: out_i = H_i[phi_i]
 *
 * MPI: Output vector gets the same MPI distribution as input vector. Only
 *      local orbitals are computed.
 */
OrbitalVector HelmholtzVector::operator()(OrbitalVector &Phi) const {
    Timer t_tot;
    Printer::printHeader(0, "Applying Helmholtz operators");
    int oldprec = Printer::setPrecision(5);

    println(0, " Orb    RealNorm   Nodes     ImagNorm   Nodes     Timing");
    Printer::printSeparator(0, '-');

    OrbitalVector out = orbital::param_copy(Phi);
    for (int i = 0; i < Phi.size(); i++) {
        if (not mpi::my_orb(out[i])) continue;

        Timer t_i;
        out[i] = apply(i, Phi[i]);
        out[i].rescale(-1.0 / (2.0 * MATHCONST::pi));
        t_i.stop();

        int rNodes = out[i].getNNodes(NUMBER::Real);
        int iNodes = out[i].getNNodes(NUMBER::Imag);
        double rNorm = 0.0;
        double iNorm = 0.0;
        if (out[i].hasReal()) rNorm = std::sqrt(out[i].real().getSquareNorm());
        if (out[i].hasImag()) iNorm = std::sqrt(out[i].imag().getSquareNorm());

        Printer::setPrecision(5);
        printout(0, std::setw(3) << i);
        printout(0, " " << std::setw(14) << rNorm);
        printout(0, " " << std::setw(5) << rNodes);
        printout(0, " " << std::setw(14) << iNorm);
        printout(0, " " << std::setw(5) << iNodes);
        printout(0, std::setw(14) << t_i.getWallTime() << std::endl);
    }

    t_tot.stop();
    Printer::printFooter(0, t_tot, 2);
    Printer::setPrecision(oldprec);
    return out;
}

OrbitalVector HelmholtzVector::operator()(const ComplexMatrix &F_mat, OrbitalVector &Phi) const {
    Timer t_tot;
    Printer::printHeader(0, "Rotating Helmholtz argument");
    ComplexMatrix L_mat = getLambdaMatrix();

    Timer rot_t;
    OrbitalVector Psi = orbital::rotate(L_mat - F_mat, Phi);
    rot_t.stop();
    Printer::printDouble(0, "Rotating orbitals", rot_t.getWallTime(), 5);

    t_tot.stop();
    Printer::printFooter(0, t_tot, 2);
    return Psi;
}

/** @brief Apply Helmholtz operator component wise on OrbitalVector
 *
 * This will construct a separate Helmholtz operator for each of the entries
 * in the OrbitalVector based on the corresponding lambda_i parameter in the
 * HelmholtzVector. Computes output as: out_i = H_i[V*phi_i + psi_i]
 *
 * Specialized version with smaller memory footprint since the full vector V*Phi
 * is never stored, but computed on the fly.
 *
 * MPI: Output vector gets the same MPI distribution as input vector. Only
 *      local orbitals are computed.
 */
OrbitalVector HelmholtzVector::operator()(RankZeroTensorOperator &V, OrbitalVector &Phi, OrbitalVector &Psi) const {
    Timer t_tot;
    Printer::printHeader(0, "Applying Helmholtz operators");
    int oldprec = Printer::setPrecision(5);

    if (Phi.size() != Psi.size()) MSG_FATAL("OrbitalVector size mismatch");

    println(0, " Orb    RealNorm   Nodes     ImagNorm   Nodes     Timing");
    Printer::printSeparator(0, '-');

    OrbitalVector out = orbital::param_copy(Phi);
    for (int i = 0; i < Phi.size(); i++) {
        if (not mpi::my_orb(out[i])) continue;

        Timer t_i;
        Orbital Vphi_i = V(Phi[i]);
        Vphi_i.add(1.0, Psi[i]);
        Vphi_i.rescale(-1.0 / (2.0 * MATHCONST::pi));
        out[i] = apply(i, Vphi_i);
        t_i.stop();

        int rNodes = out[i].getNNodes(NUMBER::Real);
        int iNodes = out[i].getNNodes(NUMBER::Imag);
        double rNorm = 0.0;
        double iNorm = 0.0;
        if (out[i].hasReal()) rNorm = std::sqrt(out[i].real().getSquareNorm());
        if (out[i].hasImag()) iNorm = std::sqrt(out[i].imag().getSquareNorm());

        Printer::setPrecision(5);
        printout(0, std::setw(3) << i);
        printout(0, " " << std::setw(14) << rNorm);
        printout(0, " " << std::setw(5) << rNodes);
        printout(0, " " << std::setw(14) << iNorm);
        printout(0, " " << std::setw(5) << iNodes);
        printout(0, std::setw(14) << t_i.getWallTime() << std::endl);
    }

    t_tot.stop();
    Printer::printFooter(0, t_tot, 2);
    Printer::setPrecision(oldprec);
    return out;
}

/** @brief Apply Helmholtz operator on individual Orbital
 *
 * This will construct a Helmholtz operator with the i-th component of the
 * lambda vector and apply it to the input orbital.
 *
 * Computes output as: out_i = H_i[phi_i]
 */
Orbital HelmholtzVector::apply(int i, Orbital &phi) const {
    ComplexDouble mu_i = std::sqrt(-2.0 * this->lambda(i));
    if (std::abs(mu_i.imag()) > mrcpp::MachineZero) MSG_FATAL("Mu cannot be complex");
    mrcpp::HelmholtzOperator H(*MRA, mu_i.real(), this->prec);

    Orbital out = phi.paramCopy();
    if (phi.hasReal()) {
        out.alloc(NUMBER::Real);
        mrcpp::apply(this->prec, out.real(), H, phi.real());
    }
    if (phi.hasImag()) {
        out.alloc(NUMBER::Imag);
        mrcpp::apply(this->prec, out.imag(), H, phi.imag());
        if (phi.conjugate()) out.imag().rescale(-1.0);
    }
    return out;
}
} // namespace mrchem
