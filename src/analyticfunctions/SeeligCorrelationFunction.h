#pragma once

#include "NuclearCorrelationFunction.h"

namespace mrchem {

class SeeligCorrelationFunction final : public NuclearCorrelationFunction {
public:
    DoubleFunction getS_m1(const Nucleus &nuc) const {
        auto f = [nuc] (const double *r) -> double {
            const double *R = nuc.getCoord();
            double r_R = math_utils::calc_distance(r, R);
            double Z = nuc.getCharge();
            return std::exp(Z*r_R);
        };
        return f;
    }
    DoubleFunction getS_0(const Nucleus &nuc) const {
        auto f = [nuc] (const double *r) -> double {
            const double *R = nuc.getCoord();
            double r_R = math_utils::calc_distance(r, R);
            double Z = nuc.getCharge();
            return std::exp(-Z*r_R);
        };
        return f;
    }
    DoubleFunction getS_1(const Nucleus &nuc, int d) const {
        auto f = [nuc, d] (const double *r) -> double {
            const double *R = nuc.getCoord();
            double r_R = math_utils::calc_distance(r, R);
            double n = (r[d]-R[d])/r_R;
            double Z = nuc.getCharge();
            return -Z*n;
        };
        return f;
    }
    DoubleFunction getS_2(const Nucleus &nuc) const {
        auto f = [nuc] (const double *r) -> double {
            double Z = nuc.getCharge();
            return -Z*Z;
        };
        return f;
    }
};

} //namespace mrchem
