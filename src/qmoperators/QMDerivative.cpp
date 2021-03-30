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

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "QMDerivative.h"

#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "utils/print_utils.h"

using mrcpp::DerivativeOperator;
using mrcpp::FunctionTree;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

QMDerivative::QMDerivative(int d, std::shared_ptr<DerivativeOperator<3>> D, bool im)
        : QMOperator()
        , imag(im)
        , apply_dir(d)
        , derivative(D) {}

Orbital QMDerivative::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");
    if (this->derivative == nullptr) MSG_ERROR("No derivative operator");

    auto dir = this->apply_dir;
    auto &D = *this->derivative;

    Orbital out = inp.paramCopy();
    if (this->isReal()) {
        if (inp.hasReal()) {
            out.alloc(NUMBER::Real);
            mrcpp::apply(out.real(), D, inp.real(), dir);
        }
        if (inp.hasImag()) {
            out.alloc(NUMBER::Imag);
            mrcpp::apply(out.imag(), D, inp.imag(), dir);
            if (inp.conjugate()) out.imag().rescale(-1.0);
        }
    } else {
        if (inp.hasImag()) {
            out.alloc(NUMBER::Real);
            mrcpp::apply(out.real(), D, inp.imag(), dir);
            if (inp.conjugate()) out.real().rescale(-1.0);
        }
        if (inp.hasReal()) {
            out.alloc(NUMBER::Imag);
            mrcpp::apply(out.imag(), D, inp.real(), dir);
            out.imag().rescale(-1.0);
        }
    }
    return out;
}

Orbital QMDerivative::dagger(Orbital inp) {
    NOT_IMPLEMENTED_ABORT;
}

} // namespace mrchem
