#pragma once

#include "FockOperator.h"

class HartreeFock : public FockOperator {
public:
    HartreeFock(KineticOperator &t,
                RegularizedPotential &v,
                CoulombOperator &j,
                ExchangeOperator &k)
            : FockOperator(&t, &v, &j, &k) { }
    virtual ~HartreeFock() { }
};

