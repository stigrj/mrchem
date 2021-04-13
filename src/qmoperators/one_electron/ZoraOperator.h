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

#pragma once

#include "tensor/RankZeroOperator.h"

#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/QMPotential.h"

namespace mrchem {

class ZoraOperator final : public RankZeroOperator {
public:
    ZoraOperator(RankZeroOperator V, double zora_factor, double proj_prec, bool mpi_share = false) {
        auto V_zora = std::make_shared<QMPotential>(1, mpi_share);
        auto f_zora = [zora_factor, &V](const mrcpp::Coord<3> &r) -> double { return zora_factor / (zora_factor - V(r).real()); };
        qmfunction::project(*V_zora, f_zora, NUMBER::Real, proj_prec);

        // Invoke operator= to assign *this operator
        RankZeroOperator &O = (*this);
        O = V_zora;
        O.name() = "V_zora";
    }
};

} // namespace mrchem
