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

#include "QMPotential.h"
#include "analyticfunctions/NuclearGradientFunction.h"
#include "qmoperators/RankOneTensorOperator.h"

namespace mrchem {

class NuclearGradientPotential final : public QMPotential {
public:
    NuclearGradientPotential(int d, double z, const mrcpp::Coord<3> &r, double c)
            : QMPotential(1)
            , func(d, z, r, c) {}

    void setup(double prec) override;
    void clear() override;

private:
    NuclearGradientFunction func;
};

class NuclearGradientOperator final : public RankOneTensorOperator<3> {
public:
    NuclearGradientOperator(double z, const mrcpp::Coord<3> &r, double c) {
        auto x_rm3 = std::make_shared<NuclearGradientPotential>(0, z, r, c);
        auto y_rm3 = std::make_shared<NuclearGradientPotential>(1, z, r, c);
        auto z_rm3 = std::make_shared<NuclearGradientPotential>(2, z, r, c);

        // Invoke operator= to assign *this operator
        RankOneTensorOperator &v = (*this);
        v[0] = x_rm3;
        v[1] = y_rm3;
        v[2] = z_rm3;
        v[0].name() = "x/r^3";
        v[1].name() = "y/r^3";
        v[2].name() = "z/r^3";
    }
};

} // namespace mrchem
