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

class DistancePotential final : public QMPotential {
public:
    DistancePotential(double pow, const mrcpp::Coord<3> &r_K, double S);

    void setup(double prec) override;
    void clear() override;

private:
    const double power;
    NuclearFunction func;
};

class DistanceOperator final : public RankZeroTensorOperator {
public:
    DistanceOperator(double pow, const mrcpp::Coord<3> &R, double S = 1.0e-7) {
        auto r_pow = std::make_shared<DistancePotential>(pow, R, S);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &v = (*this);
        v = r_pow;
        std::stringstream o_name;
        o_name << "r^{" << std::setprecision(1) << std::fixed << pow << "}";
        v.name() = o_name.str();
    }
};

} // namespace mrchem
