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

#include "qmoperators/RankOneTensorOperator.h"
#include "qmoperators/one_electron/QMPotential.h"

namespace mrchem {

class PositionPotential final : public QMPotential {
public:
    PositionPotential(int d, const mrcpp::Coord<3> &o);

protected:
    mrcpp::AnalyticFunction<3> func;

    void setup(double prec) override;
    void clear() override;
};

class PositionOperator : public RankOneTensorOperator<3> {
public:
    PositionOperator(const mrcpp::Coord<3> &o = {0.0, 0.0, 0.0}) {
        auto r_x = std::make_shared<PositionPotential>(0, o);
        auto r_y = std::make_shared<PositionPotential>(1, o);
        auto r_z = std::make_shared<PositionPotential>(2, o);

        // Invoke operator= to assign *this operator
        RankOneTensorOperator &r = (*this);
        r[0] = r_x;
        r[1] = r_y;
        r[2] = r_z;
        r[0].name() = "r[x]";
        r[1].name() = "r[y]";
        r[2].name() = "r[z]";
    }
};

} // namespace mrchem
