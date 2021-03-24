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

#include "analyticfunctions/NuclearFunction.h"
#include "qmoperators/RankZeroTensorOperator.h"
#include "qmoperators/one_electron/QMPotential.h"

namespace mrchem {

class NuclearPotential final : public QMPotential {
public:
    NuclearPotential(const Nuclei &nucs, double proj_prec, double smooth_prec = -1.0, bool mpi_share = false);
    ~NuclearPotential() override { free(NUMBER::Total); }

private:
    NuclearFunction func;

    void allreducePotential(double prec, QMFunction &V_loc);
};

class NuclearOperator final : public RankZeroTensorOperator {
public:
    NuclearOperator(const Nuclei &nucs, double proj_prec, double smooth_prec = -1.0, bool mpi_share = false) {
        auto r_m1 = std::make_shared<NuclearPotential>(nucs, proj_prec, smooth_prec, mpi_share);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &v = (*this);
        v = r_m1;
        v.name() = "V_nuc";
    }
};

} // namespace mrchem
