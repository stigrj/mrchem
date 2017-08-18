#ifndef SEELIGCORRELATIONFUNCTION_H
#define SEELIGCORRELATIONFUNCTION_H

#include "NuclearCorrelationFunction.h"

class SeeligCorrelationFunction : public NuclearCorrelationFunction {
public:
    SeeligCorrelationFunction() { }
    virtual ~SeeligCorrelationFunction() { }

    std::function<double (const double *r)> getS_m1(const Nucleus &nuc) const {
        auto f = [nuc] (const double *r) -> double {
            const double *R = nuc.getCoord();
            double r_R = MathUtils::calcDistance(3, r, R);
            double Z = nuc.getCharge();
            return 1.0/exp(-Z*r_R);
        };
        return f;
    }
    std::function<double (const double *r)> getS_0(const Nucleus &nuc) const {
        auto f = [nuc] (const double *r) -> double {
            const double *R = nuc.getCoord();
            double r_R = MathUtils::calcDistance(3, r, R);
            double Z = nuc.getCharge();
            return exp(-Z*r_R);
        };
        return f;
    }
    std::function<double (const double *r)> getS_1(const Nucleus &nuc, int d) const {
        auto f = [nuc, d] (const double *r) -> double {
            const double *R = nuc.getCoord();
            double r_R = MathUtils::calcDistance(3, r, R);
            double n = (r[d]-R[d])/r_R;
            double Z = nuc.getCharge();
            return -Z*n;
        };
        return f;
    }
    std::function<double (const double *r)> getS_2(const Nucleus &nuc) const {
        auto f = [nuc] (const double *r) -> double {
            double Z = nuc.getCharge();
            return -Z*Z;
        };
        return f;
    }
};

#endif // SEELIGCORRELATIONFUNCTION_H
