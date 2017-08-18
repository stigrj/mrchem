#ifndef IDENTITYCORRELATIONFUNCTION_H
#define IDENTITYCORRELATIONFUNCTION_H

#include "NuclearCorrelationFunction.h"

class IdentityCorrelationFunction : public NuclearCorrelationFunction {
public:
    IdentityCorrelationFunction() { }
    virtual ~IdentityCorrelationFunction() { }

    std::function<double (const double *r)> getS_m1(const Nucleus &nuc) const {
        auto f = [nuc] (const double *r) -> double { return 1.0; };
        return f;
    }
    std::function<double (const double *r)> getS_0(const Nucleus &nuc) const {
        auto f = [nuc] (const double *r) -> double { return 1.0; };
        return f;
    }
    std::function<double (const double *r)> getS_1(const Nucleus &nuc, int d) const {
        auto f = [nuc, d] (const double *r) -> double { return 0.0; };
        return f;
    }
    std::function<double (const double *r)> getS_2(const Nucleus &nuc) const {
        auto f = [nuc] (const double *r) -> double {
            const double *R = nuc.getCoord();
            double r_R = MathUtils::calcDistance(3, r, R);
            double Z = nuc.getCharge();
            return Z/r_R;
        };
        return f;
    }
};

#endif // IDENTITYCORRELATIONFUNCTION_H
