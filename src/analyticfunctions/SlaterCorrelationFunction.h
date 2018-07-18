#pragma once

#include "NuclearCorrelationFunction.h"

namespace mrchem {

class SlaterCorrelationFunction final : public NuclearCorrelationFunction {
public:
    SlaterCorrelationFunction(double a) : param(a) { }

    DoubleFunction getS_m1(const Nucleus &nuc) const {
        double a = this->param;
        auto f = [a, nuc] (const double *r) -> double {
            const double *R = nuc.getCoord();
            double r_R = math_utils::calc_distance(r, R);
            double Z = nuc.getCharge();
            return 1.0/(1.0 + std::exp(-a*r_R*Z)/(a - 1.0));
        };
        return f;
    }
    DoubleFunction getS_0(const Nucleus &nuc) const {
        double a = this->param;
        auto f = [a, nuc] (const double *r) -> double {
            const double *R = nuc.getCoord();
            double r_R = math_utils::calc_distance(r, R);
            double Z = nuc.getCharge();
            return 1.0 + std::exp(-a*r_R*Z)/(a - 1.0);
        };
        return f;
    }
    DoubleFunction getS_1(const Nucleus &nuc, int d) const {
        double a = this->param;
        auto f = [a, nuc, d] (const double *r) -> double {
            const double *R = nuc.getCoord();
            double r_R = math_utils::calc_distance(r, R);
            double n = (r[d]-R[d])/r_R;
            double Z = nuc.getCharge();
            return -(a*Z)/(1.0 + (a - 1.0)*std::exp(a*r_R*Z))*n;
        };
        return f;
    }
    DoubleFunction getS_2(const Nucleus &nuc) const {
        double a = this->param;
        auto f = [a, nuc] (const double *r) -> double {
            const double *R = nuc.getCoord();
            double r_R = math_utils::calc_distance(r, R);
            double Z = nuc.getCharge();
            double nomin = -1.0;
            double denom = r_R;
            if (r_R < 16.0) {
                nomin = 2.0*(a - 1.0)*(1.0 - std::exp(a*r_R*Z)) - a*a*r_R*Z;
                denom = 2.0*r_R*(1.0 + (a - 1.0)*std::exp(a*r_R*Z));
            }
            return -Z*nomin/denom;
        };
        return f;
    }
protected:
    const double param;
};

} //namespace mrchem
