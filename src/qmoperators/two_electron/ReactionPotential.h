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

#include "qmoperators/QMPotential.h"

#include "environment/SCRF.h"

namespace mrchem {
/** @class ReactionPotential
 *  @brief class containing the solvent-substrate interaction reaction potential
 *  obtained by solving
 *  \f[
 *     \Delta V_{R} = -4\pi\left( \rho\frac{1-\epsilon}{\epsilon} + \gamma_s  \right)
 *  \f]
 *  where \f$\rho\f$ is the total molecular density of a solute molecule, \f$\epsilon\f$ is
 *  the Permittivity function of the continuum and \f$\gamma_s\f$ is the surface charge distribution.
 */
class ReactionPotential final : public QMPotential {
public:
    /** @brief Initializes the ReactionPotential class.
     *  @param scrf_p A SCRF instance which contains the parameters needed to compute the ReactionPotential.
     *  @param Phi_p A pointer to a vector which contains the orbitals optimized in the SCF procedure. */
    ReactionPotential(std::unique_ptr<SCRF> &scrf_p, std::shared_ptr<mrchem::OrbitalVector> Phi_p);
    ~ReactionPotential() override { free(NUMBER::Total); }

    SCRF *getHelper() { return this->helper.get(); }
    double getElectronicEnergy() { return this->helper->getElectronicEnergy(); }
    double getNuclearEnergy() { return this->helper->getNuclearEnergy(); }
    double getTotalEnergy() { return this->helper->getTotalEnergy(); }

    /** @brief Updates the helper.mo_residual member variable. This variable is used to set the convergence criterion in
     * the dynamic convergence method. */
    void updateMOResidual(double const err_t) { this->helper->mo_residual = err_t; }

    QMFunction &getCurrentReactionPotential() { return this->helper->getCurrentReactionPotential(); }
    QMFunction &getPreviousReactionPotential() { return this->helper->getPreviousReactionPotential(); }
    QMFunction &getCurrentDifferenceReactionPotential() { return this->helper->getCurrentDifferenceReactionPotential(); }

    QMFunction &getCurrentGamma() { return this->helper->getCurrentGamma(); }
    QMFunction &getPreviousGamma() { return this->helper->getPreviousGamma(); }
    QMFunction &getCurrentDifferenceGamma() { return this->helper->getCurrentDifferenceGamma(); }

    /** @brief changes the #first_iteration variable to false. This is done in tests where we want the ReactionPotential
     * to be calculated at the zero-th iteration in the SCF cycle instead of not doing it, as is the deafult.*/
    void setTesting() { this->first_iteration = false; }

    friend class ReactionOperator;

protected:
    void clear();

private:
    bool first_iteration = true; //!< when this variable is true the reaction potential is not calculated. This is done
                                 //!< only in the zero-th iteration of the SCF procedure, after that it is set to false.
    std::unique_ptr<SCRF> helper; //!< A SCRF instance used to compute the ReactionPotential.
    std::shared_ptr<mrchem::OrbitalVector> Phi; //!< holds the Orbitals needed to compute the electronic density for the SCRF procedure.

    void setup(double prec);
};

} // namespace mrchem
