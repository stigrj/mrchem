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

#include "ZoraPotential.h"
#include "analyticfunctions/NuclearFunction.h"
#include "qmoperators/QMPotential.h"
#include "tensor/RankZeroOperator.h"

using DerivativeOperator_p = std::shared_ptr<mrcpp::DerivativeOperator<3>>;

namespace mrchem {

class ZoraOperator final : public RankZeroOperator {
public:
    ZoraOperator(const Nuclei &nucs, double zora_factor, double proj_prec, double smooth_prec = -1.0, bool mpi_share = false) {
        auto kappa = std::make_shared<ZoraPotential>(nucs, zora_factor, proj_prec, smooth_prec, mpi_share);

        // Invoke operator= to assign *this operator
        RankZeroOperator &v = (*this);
        v = kappa;
        v.name() = "V_zora";
    }
};

} // namespace mrchem
