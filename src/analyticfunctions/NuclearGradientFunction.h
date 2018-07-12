#pragma once

#include "MRCPP/MWFunctions"

#include "Nucleus.h"

namespace mrchem {

class NuclearGradientFunction : public mrcpp::RepresentableFunction<3> {
public:
    NuclearGradientFunction(int d, const Nucleus &nuc, double c)
        : dir(d), smooth(c), nucleus(nuc) { }

    double evalf(const double *r) const;

    Nucleus &getNucleus() { return this->nucleus; }
    const Nucleus &getNucleus() const { return this->nucleus; }

    bool isVisibleAtScale(int scale, int nQuadPts) const;
    bool isZeroOnInterval(const double *a, const double *b) const;

protected:
    int dir;
    double smooth;
    Nucleus nucleus;

    double du_dr(double r1) const;
};

} //namespace mrchem
