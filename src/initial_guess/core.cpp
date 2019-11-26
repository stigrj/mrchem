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

#include <Eigen/Eigenvalues>

#include "MRCPP/MWFunctions"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"
#include "utils/math_utils.h"

#include "core.h"

#include "analyticfunctions/HydrogenFunction.h"
#include "chemistry/Molecule.h"
#include "chemistry/Nucleus.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/OrbitalIterator.h"
#include "qmfunctions/orbital_utils.h"
#include "qmfunctions/qmfunction_utils.h"

#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

namespace initial_guess {
namespace core {
/** @brief Helper struct to get the orbital ordering right
 *
 *  First index energy level (n)
 *  Second index angular momentum (l)
 */
// clang-format off
int PT[29][2] = {
     /*s*/
    {1, 0},                      /*p*/
    {2, 0},                     {2, 1},
    {3, 0},               /*d*/ {3, 1},
    {4, 0},              {3, 2},{4, 1},
    {5, 0},        /*f*/ {4, 2},{5, 1},
    {6, 0},       {4, 3},{5, 2},{6, 1},
    {7, 0}, /*g*/ {5, 3},{6, 2},{7, 1},
    {8, 0},{5, 4},{6, 3},{7, 2},{8, 1},
    {9, 0},{6, 4},{7, 3},{8, 2},{9, 1}
};
// clang-format on

} // namespace core
} // namespace initial_guess

/** @brief Produce an initial guess of orbitals
 *
 * @param prec: precision used in projection
 * @param mol: molecule
 * @param restricted: spin restriction
 * @param zeta: quality of hydrogen AO basis
 *
 * Sets up an AO basis of hydrogen functions with the given zeta quality
 * (SZ, DZ, TZ, QZ), computes and diagonalizes the core Hamiltonian matrix,
 * and fills the resulting orbitals by the Aufbau principle.
 *
 */
OrbitalVector initial_guess::core::setup(double prec, const Molecule &mol, bool restricted, int zeta) {
    std::stringstream o_prec, o_zeta;
    o_prec << std::setprecision(5) << std::scientific << prec;
    o_zeta << zeta;
    mrcpp::print::separator(0, '~');
    print_utils::text(0, "Calculation ", "Diagonalize Hamiltonian matrix");
    print_utils::text(0, "Precision   ", o_prec.str());
    print_utils::text(0, "Restricted  ", (restricted) ? "True" : "False");
    print_utils::text(0, "Hamiltonian ", "Core");
    print_utils::text(0, "AO basis    ", "Hydrogenic orbitals");
    print_utils::text(0, "Zeta quality", o_zeta.str());
    mrcpp::print::separator(0, '~', 2);

    Timer t_tot, t_lap;
    auto plevel = Printer::getPrintLevel();
    if (plevel == 1) mrcpp::print::header(1, "Core-Hamiltonian Initial Guess");

    int mult = mol.getMultiplicity(); // multiplicity
    int Ne = mol.getNElectrons();     // total electrons
    int Nd = Ne - (mult - 1);         // doubly occupied
    if (Nd % 2 != 0) MSG_ABORT("Invalid multiplicity");

    // Make Fock operator contributions
    t_lap.start();
    auto D_p = std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5);
    KineticOperator T(D_p);
    NuclearOperator V(mol.getNuclei(), prec);
    if (plevel == 1) mrcpp::print::time(1, "Projecting nuclear potential", t_lap);

    // Project AO basis of hydrogen functions
    t_lap.start();
    OrbitalVector Phi = initial_guess::core::project_ao(prec, mol.getNuclei(), SPIN::Paired, zeta);
    ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(Phi);
    if (plevel == 1) mrcpp::print::time(1, "Projecting Hydrogen AOs", t_lap);

    // Compute Hamiltonian matrix
    Timer t_diag;
    mrcpp::print::header(2, "Diagonalize Core-Hamiltonian matrix");

    t_lap.start();
    T.setup(prec);
    V.setup(prec);
    mrcpp::print::time(1, "Building Fock operator", t_lap);
    t_lap.start();
    ComplexMatrix t = T(Phi, Phi);
    ComplexMatrix v = V(Phi, Phi);
    ComplexMatrix F = S_m12.transpose() * (t + v) * S_m12;
    V.clear();
    T.clear();
    mrcpp::print::time(1, "Computing Fock matrix", t_lap);

    // Diagonalize Hamiltonian matrix
    t_lap.start();
    Eigen::SelfAdjointEigenSolver<ComplexMatrix> es(F.cols());
    es.compute(F);
    ComplexMatrix ei_vec = es.eigenvectors();
    ComplexMatrix U = ei_vec.transpose() * S_m12;
    mrcpp::print::time(1, "Diagonalizing Fock matrix", t_lap);

    // Need to convert to QMFunctions for linear_combination
    QMFunctionVector funcs;
    for (auto &phi_i : Phi) funcs.push_back(phi_i);

    // Rotate orbitals and fill electrons by Aufbau
    t_lap.start();
    OrbitalVector Psi;
    if (restricted) {
        if (mult != 1) MSG_ABORT("Restricted open-shell not available");
        int Np = Nd / 2; // paired orbitals
        Psi = initial_guess::core::rotate_orbitals(prec, U, Phi, Np, SPIN::Paired);
    } else {
        int Na = Nd / 2 + (mult - 1); // alpha orbitals
        int Nb = Nd / 2;              // beta orbitals
        OrbitalVector Psi_a = initial_guess::core::rotate_orbitals(prec, U, Phi, Na, SPIN::Alpha);
        OrbitalVector Psi_b = initial_guess::core::rotate_orbitals(prec, U, Phi, Nb, SPIN::Beta);
        Psi = orbital::adjoin(Psi_a, Psi_b);
    }
    mrcpp::print::time(1, "Rotating orbitals", t_lap);
    mrcpp::print::footer(2, t_diag, 2);
    if (plevel == 1) mrcpp::print::footer(1, t_tot, 2);

    return Psi;
}

/** @brief Project AO basis of hydrogen functions
 *
 * @param prec: precision used in projection
 * @param nucs: the nuclei of the molecule
 * @param zeta: quality of hydrogen AO basis
 *
 * Sets up an AO basis of hydrogen functions with the given zeta quality
 * (SZ, DZ, TZ, QZ), and projects it onto the MW basis. The basis at each
 * atomic center is always a complete shell: single zeta (SZ) means that
 * the current shell is completed, double zeta (DZ) means that also the
 * next shell will be completed, etc. E.i. the oxygen atom will get the
 * following AO basis:
 *
 * Oxygen AOs:
 * SZ: 1s2s2p                 (2s +  3p)
 * DZ: 1s2s2p3s3p             (3s +  6p)
 * TZ: 1s2s2p3s3p4s3d4p       (4s +  9p +  5d)
 * QZ: 1s2s2p3s3p4s3d4p5s4d5p (5s + 12p + 10d)
 *
 */
OrbitalVector initial_guess::core::project_ao(double prec, const Nuclei &nucs, int spin, int zeta) {
    Timer t_tot;
    auto w0 = Printer::getWidth() - 2;
    auto w1 = 5;
    auto w2 = 7;
    auto w3 = w0 * 2 / 9;
    auto w4 = w0 - w1 - w2 - 3 * w3;

    std::stringstream o_head;
    o_head << std::setw(w1) << "n";
    o_head << std::setw(w4) << "Atom";
    o_head << std::setw(w2) << "Label";
    o_head << std::setw(w3 + 1) << "Nodes";
    o_head << std::setw(w3) << "Size";
    o_head << std::setw(w3) << "Time";

    mrcpp::print::header(2, "Projecting Hydrogen AOs");
    println(2, o_head.str());
    mrcpp::print::separator(2, '-');

    const char label[10] = "spdfg";

    OrbitalVector Phi;
    for (int i = 0; i < nucs.size(); i++) {
        const Nucleus &nuc = nucs[i];
        int minAO = std::ceil(nuc.getElement().getZ() / 2.0);
        double Z = nuc.getCharge();
        const mrcpp::Coord<3> &R = nuc.getCoord();

        int nAO = 0;
        int nShell = 0;
        int zetaReached = 0;
        bool minAOReached = false;
        while (true) {
            int n = initial_guess::core::PT[nShell][0];
            int l = initial_guess::core::PT[nShell][1];
            int M = 2 * l + 1;

            if (minAOReached and l == 0) zetaReached++;
            if (zetaReached >= zeta) break;

            for (int m = 0; m < M; m++) {
                Timer t_i;
                HydrogenFunction h_func(n, l, m, Z, R);

                Phi.push_back(Orbital(spin));
                Phi.back().setRankID(Phi.size() % mpi::orb_size);
                if (mpi::my_orb(Phi.back())) qmfunction::project(Phi.back(), h_func, NUMBER::Real, prec);

                std::stringstream o_txt;
                o_txt << std::setw(w1 - 1) << Phi.size() - 1;
                o_txt << std::setw(w4) << nuc.getElement().getSymbol();
                o_txt << std::setw(w2 - 1) << n << label[l];
                print_utils::qmfunction(2, o_txt.str(), Phi.back(), t_i);

                if (++nAO >= minAO) minAOReached = true;
            }
            nShell++;
        }
    }
    mrcpp::print::footer(2, t_tot, 2);
    return Phi;
}

OrbitalVector initial_guess::core::rotate_orbitals(double prec, ComplexMatrix &U, OrbitalVector &Phi, int N, int spin) {
    Timer t_tot;
    OrbitalVector Psi;
    for (int i = 0; i < N; i++) Psi.push_back(Orbital(spin));
    mpi::distribute(Psi);

    OrbitalIterator iter(Phi);
    while (iter.next()) {
        for (int i = 0; i < Psi.size(); i++) {
            if (not mpi::my_orb(Psi[i])) continue;
            QMFunctionVector func_vec;
            ComplexVector coef_vec(iter.get_size());
            for (int j = 0; j < iter.get_size(); j++) {
                int idx_j = iter.idx(j);
                Orbital &recv_j = iter.orbital(j);
                coef_vec[j] = U(i, idx_j);
                func_vec.push_back(recv_j);
            }
            Orbital tmp_i = Psi[i].paramCopy();
            qmfunction::linear_combination(tmp_i, coef_vec, func_vec, prec);
            Psi[i].add(1.0, tmp_i); // In place addition
            Psi[i].crop(prec);
        }
    }
    mrcpp::print::time(1, "Rotating orbitals", t_tot);
    return Psi;
}

} // namespace mrchem
