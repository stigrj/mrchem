#pragma once

#include "CoulombPotential.h"

namespace mrchem {

class CoulombPotentialD2 final : public CoulombPotential {
public:
    CoulombPotentialD2(mrcpp::PoissonOperator *P, OrbitalVector *Phi, OrbitalVector *X, OrbitalVector *Y);

private:
    OrbitalVector *orbitals_x; ///< Perturbed orbitals
    OrbitalVector *orbitals_y; ///< Perturbed orbitals

    void setupDensity(double prec);
};

} //namespace mrchem
