#pragma once

#include "MWAdder.h"
#include "MWMultiplier.h"
#include "GridGenerator.h"

class OrbitalVector;
class Orbital;
class Density;
class QMOperator;

class DensityProjector {
public:
    DensityProjector(double prec, int max_scale);
    virtual ~DensityProjector() { }

    void setPrecision(double prec);

    void operator()(Density &rho, Orbital &phi, QMOperator *R = 0);
    void operator()(Density &rho, OrbitalVector &phi, QMOperator *R = 0);

protected:
    MWAdder<3> add;
    MWMultiplier<3> mult;
    GridGenerator<3> grid;
};

