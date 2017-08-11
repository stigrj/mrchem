#pragma once

#include "QMPotential.h"
#include "Density.h"

class CoulombOperator : public QMPotential {
public:
    CoulombOperator(PoissonOperator &P, OrbitalVector &phi, QMOperator *R)
        : poisson(&P),
          orbitals(&phi),
          nuc_corr_fac(R),
#ifdef HAVE_MPI
	density(false, true){//Use shared memory.
	//          density(false, false){//do not Use shared memory.
#else
          density(false) {
#endif
    }
    virtual ~CoulombOperator() { }

protected:
    PoissonOperator *poisson;   // Pointer to external object
    OrbitalVector *orbitals;    // Pointer to external object
    QMOperator *nuc_corr_fac;   // Nuclear correlation factor
    Density density;            // Density that defines the potential
};

