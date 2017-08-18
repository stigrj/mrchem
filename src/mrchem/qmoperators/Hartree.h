#pragma once

#include "FockOperator.h"

class Hartree : public FockOperator {
public:
    Hartree(KineticOperator &t,
            RegularizedPotential &v,
            CoulombOperator &j)
        : FockOperator(&t, &v, &j) { }
    virtual ~Hartree() { }
};

