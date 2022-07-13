/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "Permittivity.h"
#include "Cavity.h"
#include <MRCPP/MWFunctions>

namespace mrchem {

Permittivity::Permittivity(const mrchem::Cavity cavity, double epsilon_in, double epsilon_out, std::string formulation)
        : epsilon_in(epsilon_in)
        , epsilon_out(epsilon_out)
        , formulation(formulation)
        , cavity(cavity) {}

double Permittivity::evalf(const mrcpp::Coord<3> &r) const {
    auto epsilon = epsilon_in * std::exp(std::log(epsilon_out / epsilon_in) * (1 - this->cavity.evalf(r)));
    if (inverse) {
        return 1 / epsilon;
    } else {
        return epsilon;
    }
}

void Permittivity::printParameters() {
    // Collect relevant quantities
    Cavity cavity = getCavity();
    std::vector<mrcpp::Coord<3>> coords = cavity.getCoordinates();
    std::vector<double> radii = cavity.getRadii();

    // Set widths
    auto w0 = mrcpp::Printer::getWidth() - 1;
    auto w1 = 5;
    auto w2 = 12;
    auto w3 = 2 * w0 / 9;
    auto w4 = w0 - w1 - w2 - 3 * w3;

    // Build table column headers
    std::stringstream o_head;
    o_head << std::setw(w1) << "N";
    o_head << std::setw(w2) << "Radius";
    o_head << std::string(w4 - 1, ' ') << ':';
    o_head << std::setw(w3) << "x";
    o_head << std::setw(w3) << "y";
    o_head << std::setw(w3) << "z";

    // Print
    mrcpp::print::header(0, "Solvation Cavity");
    print_utils::text(0, "Formulation", getFormulation(), true);
    print_utils::scalar(0, "Cavity width", cavity.getWidth(), "", 6);
    print_utils::scalar(0, "Dielectric constant", getEpsIn(), "(in)", 6);
    print_utils::scalar(0, "", getEpsOut(), "(out)", 6);
    mrcpp::print::separator(0, '-');
    println(0, o_head.str());
    mrcpp::print::separator(0, '-');
    for (int i = 0; i < coords.size(); i++) {
        mrcpp::Coord<3> coord = coords[i];
        double r = radii[i];
        double x = coord[0];
        double y = coord[1];
        double z = coord[2];

        std::stringstream o_coord;
        o_coord << std::setw(w1) << i;
        o_coord << std::setw(w2) << std::setprecision(6) << std::fixed << r;
        o_coord << std::string(w4 - 1, ' ') << ':';
        o_coord << std::setw(w3) << std::setprecision(6) << std::fixed << x;
        o_coord << std::setw(w3) << std::setprecision(6) << std::fixed << y;
        o_coord << std::setw(w3) << std::setprecision(6) << std::fixed << z;
        println(0, o_coord.str());
    }
    mrcpp::print::separator(0, '=', 2);
}

} // namespace mrchem
