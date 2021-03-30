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

#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "QMPotential.h"

#include "qmfunctions/Orbital.h"
#include "utils/print_utils.h"

using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief constructor
 *
 * @param adap: extra refinement in output
 *
 * Initializes the QMFunction with NULL pointers for both real and imaginary part.
 * These must be computed in setup() of derived classes. The initial output grid
 * in application will be a copy of the input orbital but NOT a copy of the
 * potential grid. The argument sets how many extra refinement levels is allowed
 * beyond this initial refinement.
 */
QMPotential::QMPotential(int adap, bool shared)
        : QMFunction(shared)
        , QMOperator()
        , adap_build(adap) {}

/** @brief apply potential
 *
 * @param inp: orbital on which to apply
 *
 * Computes a new orbital that is the product of this potential and the input
 * orbital. Orbital parameters are copied from input.
 */
Orbital QMPotential::apply(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    Orbital out = inp.paramCopy();
    calcRealPart(out, inp, false);
    calcImagPart(out, inp, false);

    return out;
}

/** @brief apply complex cojugate potential
 *
 * @param inp: orbital on which to apply
 *
 * Computes a new orbital that is the product of the complex conjugate of this
 * potential and the input orbital. Orbital parameters are copied from input.
 */
Orbital QMPotential::dagger(Orbital inp) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    Orbital out = inp.paramCopy();
    calcRealPart(out, inp, true);
    calcImagPart(out, inp, true);

    return out;
}

/** @brief compute real part of output
 *
 * @param inp: input orbital
 * @param dagger: apply complex conjugate potential
 *
 * Computes the real part of the output orbital. The initial output grid is a
 * copy of the input orbital grid but NOT a copy of the potential grid.
 */
void QMPotential::calcRealPart(Orbital &out, Orbital &inp, bool dagger) {
    int adap = this->adap_build;
    double prec = this->apply_prec;

    if (out.hasReal()) MSG_ABORT("Output not empty");
    if (out.isShared()) MSG_ABORT("Cannot share this function");

    QMFunction &V = *this;
    if (V.hasReal() and inp.hasReal()) {
        double coef = 1.0;
        Orbital tmp = out.paramCopy();
        tmp.alloc(NUMBER::Real);
        mrcpp::copy_grid(tmp.real(), inp.real());
        mrcpp::multiply(prec, tmp.real(), coef, V.real(), inp.real(), adap);
        out.add(1.0, tmp);
    }
    if (V.hasImag() and inp.hasImag()) {
        double coef = -1.0;
        if (dagger) coef *= -1.0;
        if (inp.conjugate()) coef *= -1.0;
        Orbital tmp = out.paramCopy();
        tmp.alloc(NUMBER::Real);
        mrcpp::copy_grid(tmp.real(), inp.imag());
        mrcpp::multiply(prec, tmp.real(), coef, V.imag(), inp.imag(), adap);
        out.add(1.0, tmp);
    }
}

/** @brief compute imaginary part of output
 *
 * @param inp: input orbital
 * @param dagger: apply complex conjugate potential
 *
 * Computes the imaginary part of the output orbital. The initial output grid is a
 * copy of the input orbital grid but NOT a copy of the potential grid.
 */
void QMPotential::calcImagPart(Orbital &out, Orbital &inp, bool dagger) {
    int adap = this->adap_build;
    double prec = this->apply_prec;

    if (out.hasImag()) MSG_ABORT("Output not empty");
    if (out.isShared()) MSG_ABORT("Cannot share this function");

    QMFunction &V = *this;
    if (V.hasReal() and inp.hasImag()) {
        double coef = 1.0;
        if (inp.conjugate()) coef *= -1.0;
        Orbital tmp = out.paramCopy();
        tmp.alloc(NUMBER::Imag);
        mrcpp::copy_grid(tmp.imag(), inp.imag());
        mrcpp::multiply(prec, tmp.imag(), coef, V.real(), inp.imag(), adap);
        out.add(1.0, tmp);
    }
    if (V.hasImag() and inp.hasReal()) {
        double coef = 1.0;
        if (dagger) coef *= -1.0;
        Orbital tmp = out.paramCopy();
        tmp.alloc(NUMBER::Imag);
        mrcpp::copy_grid(tmp.imag(), inp.real());
        mrcpp::multiply(prec, tmp.imag(), coef, V.imag(), inp.real(), adap);
        out.add(1.0, tmp);
    }
}

} // namespace mrchem
