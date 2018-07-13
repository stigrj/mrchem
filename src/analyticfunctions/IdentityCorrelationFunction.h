#pragma once

#include "NuclearCorrelationFunction.h"

namespace mrchem {

class IdentityCorrelationFunction final : public NuclearCorrelationFunction {
public:
    DoubleFunction getS_m1(const Nucleus &nuc) const {
        auto f = [nuc] (const double *r) -> double { return 1.0; };
        return f;
    }
    DoubleFunction getS_0(const Nucleus &nuc) const {
        auto f = [nuc] (const double *r) -> double { return 1.0; };
        return f;
    }
    DoubleFunction getS_1(const Nucleus &nuc, int d) const {
        auto f = [nuc, d] (const double *r) -> double { return 0.0; };
        return f;
    }
    DoubleFunction getS_2(const Nucleus &nuc) const {
        auto f = [nuc] (const double *r) -> double {
            const double *R = nuc.getCoord();
            double r_R = math_utils::calc_distance(r, R);
            double Z = nuc.getCharge();
            return Z/r_R;
        };
        return f;
    }
};

} //namespace mrchem
