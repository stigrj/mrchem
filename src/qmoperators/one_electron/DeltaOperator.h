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
#include "qmoperators/RankZeroTensorOperator.h"

namespace mrchem {

class QMDelta final : public QMPotential {
public:
    QMDelta(const mrcpp::Coord<3> &o, double expo);

private:
    mrcpp::GaussFunc<3> func;

    void setup(double prec) override;
    void clear() override;
};

class DeltaOperator final : public RankZeroTensorOperator {
public:
    DeltaOperator(const mrcpp::Coord<3> &o, double expo = 1.0e6) {
        auto delta = std::make_shared<QMDelta>(o, expo);

        // Invoke operator= to assign *this operator
        RankZeroTensorOperator &d = (*this);
        d = delta;
        d.name() = "delta";
    }
};

} // namespace mrchem
