#pragma once

#include "CoulombOperator.h"

class CoulombPotential : public CoulombOperator {
public:
    CoulombPotential(PoissonOperator &P, OrbitalVector &phi, QMOperator *R = 0)
        : CoulombOperator(P, phi, R) { }
    virtual ~CoulombPotential() { }

    virtual void setup(double prec) {
        setApplyPrec(prec);
        calcDensity(this->density, *this->orbitals, this->nuc_corr_fac);
        calcPotential(*this, this->density);
    }
    virtual void clear() {
        clearReal(true);
        clearImag(true);
        this->density.clear();
        clearApplyPrec();
    }

protected:
    void calcDensity(Density &rho, OrbitalVector &phi, QMOperator *R);
    void calcPotential(QMPotential &V, Density &rho);
};

