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

#include "ReactionPotential.h"

#include "qmfunctions/qmfunction_utils.h"

using SCRF_p = std::unique_ptr<mrchem::SCRF>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {

ReactionPotential::ReactionPotential(SCRF_p scrf_p, OrbitalVector_p Phi_p)
        : QMPotential(1, false)
        , helper(std::move(scrf_p))
        , Phi(Phi_p) {}

void ReactionPotential::setup(double prec) {
    setApplyPrec(prec);
    if (this->first_iteration) {
        this->first_iteration = false;
        return;
    }
    auto potential = this->helper->setup(prec, this->Phi);
    qmfunction::deep_copy(*this, potential);
}

void ReactionPotential::clear() {
    clearApplyPrec();
    this->helper->clear();
}

} // namespace mrchem
