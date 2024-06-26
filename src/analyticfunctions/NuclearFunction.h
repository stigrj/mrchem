/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2023 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include <cmath>
#include <vector>

#include <MRCPP/MWFunctions>

#include "chemistry/Nucleus.h"

namespace mrchem {

class NuclearFunction : public mrcpp::RepresentableFunction<3> {
public:
    NuclearFunction() = default;
    void push_back(const std::string &atom, const mrcpp::Coord<3> &r, double p1, double p2);
    void push_back(const Nucleus &nuc, double p1, double p2);

    double getPrec() const { return this->prec; }
    Nuclei &getNuclei() { return this->nuclei; }
    const Nuclei &getNuclei() const { return this->nuclei; }

    bool isVisibleAtScale(int scale, int nQuadPts) const override;
    bool isZeroOnInterval(const double *a, const double *b) const override;

    virtual std::string getParamName1() const = 0;
    virtual std::string getParamName2() const = 0;
    virtual double calcParam1(double prec, const Nucleus &nuc) const = 0;
    virtual double calcParam2(double prec, const Nucleus &nuc) const = 0;

protected:
    double prec;
    Nuclei nuclei;
    std::vector<double> param1;
    std::vector<double> param2;
};

} // namespace mrchem
