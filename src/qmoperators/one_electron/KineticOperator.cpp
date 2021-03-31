/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2021 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "KineticOperator.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief Expectation value matrix
 *
 * @param bra: orbitals on the lhs
 * @param ket: orbitals on the rhs
 *
 * Instead of applying the full kinetic operator on the ket's, the momentum
 * operator is applied both to the left and right, thus taking advantage
 * of symmetry and getting away with only first-derivative operators.
 */
ComplexMatrix qmoperator::calc_kinetic_matrix(KineticOperator &T, OrbitalVector &bra, OrbitalVector &ket) {
    RankZeroOperator p_x = T.get(0, 0);
    RankZeroOperator p_y = T.get(1, 0);
    RankZeroOperator p_z = T.get(2, 0);

    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix T_x = ComplexMatrix::Zero(Ni, Nj);
    ComplexMatrix T_y = ComplexMatrix::Zero(Ni, Nj);
    ComplexMatrix T_z = ComplexMatrix::Zero(Ni, Nj);
    {
        Timer timer;
        int nNodes = 0, sNodes = 0;
        if (&bra == &ket) {
            OrbitalVector dKet = p_x(ket);
            nNodes += orbital::get_n_nodes(dKet);
            sNodes += orbital::get_size_nodes(dKet);
            T_x = orbital::calc_overlap_matrix(dKet);
        } else {
            OrbitalVector dBra = p_x(bra);
            OrbitalVector dKet = p_x(ket);
            nNodes += orbital::get_n_nodes(dBra);
            nNodes += orbital::get_n_nodes(dKet);
            sNodes += orbital::get_size_nodes(dBra);
            sNodes += orbital::get_size_nodes(dKet);
            T_x = orbital::calc_overlap_matrix(dBra, dKet);
        }
        mrcpp::print::tree(2, "<i|p[x]p[x]|j>", nNodes, sNodes, timer.elapsed());
    }
    {
        Timer timer;
        int nNodes = 0, sNodes = 0;
        if (&bra == &ket) {
            OrbitalVector dKet = p_y(ket);
            nNodes += orbital::get_n_nodes(dKet);
            sNodes += orbital::get_size_nodes(dKet);
            T_y = orbital::calc_overlap_matrix(dKet);
        } else {
            OrbitalVector dBra = p_y(bra);
            OrbitalVector dKet = p_y(ket);
            nNodes += orbital::get_n_nodes(dBra);
            nNodes += orbital::get_n_nodes(dKet);
            sNodes += orbital::get_size_nodes(dBra);
            sNodes += orbital::get_size_nodes(dKet);
            T_y = orbital::calc_overlap_matrix(dBra, dKet);
        }
        mrcpp::print::tree(2, "<i|p[y]p[y]|j>", nNodes, sNodes, timer.elapsed());
    }
    {
        Timer timer;
        int nNodes = 0, sNodes = 0;
        if (&bra == &ket) {
            OrbitalVector dKet = p_z(ket);
            nNodes += orbital::get_n_nodes(dKet);
            sNodes += orbital::get_size_nodes(dKet);
            T_z = orbital::calc_overlap_matrix(dKet);
        } else {
            OrbitalVector dBra = p_z(bra);
            OrbitalVector dKet = p_z(ket);
            nNodes += orbital::get_n_nodes(dBra);
            nNodes += orbital::get_n_nodes(dKet);
            sNodes += orbital::get_size_nodes(dBra);
            sNodes += orbital::get_size_nodes(dKet);
            T_z = orbital::calc_overlap_matrix(dBra, dKet);
        }
        mrcpp::print::tree(2, "<i|p[z]p[z]|j>", nNodes, sNodes, timer.elapsed());
    }

    return 0.5 * (T_x + T_y + T_z);
}

} // namespace mrchem
