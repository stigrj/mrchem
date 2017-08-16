#ifndef SLATERCORRELATIONFUNCTION_H
#define SLATERCORRELATIONFUNCTION_H

#include "NuclearCorrelationFunction.h"

class SlaterCorrelationFunction : public NuclearCorrelationFunction {
public:
    SlaterCorrelationFunction(double a) : param(a) { }
    virtual ~SlaterCorrelationFunction() { }

    std::function<double (const double *r)> getS_0(const Nucleus &nuc) const {
        double a = this->param;
        auto f = [a, nuc] (const double *r) -> double {
            const double *R = nuc.getCoord();
            double r_R = MathUtils::calcDistance(3, r, R);
            double Z = nuc.getCharge();
            return 1.0 + exp(-a*r_R*Z)/(a - 1.0);
        };
        return f;
    }
    std::function<double (const double *r)> getS_1(const Nucleus &nuc, int d) const {
        double a = this->param;
        auto f = [a, nuc, d] (const double *r) -> double {
            const double *R = nuc.getCoord();
            double r_R = MathUtils::calcDistance(3, r, R);
            double n = (r[d]-R[d])/r_R;
            double Z = nuc.getCharge();
            return -(a*Z)/(1.0 + (a - 1.0)*exp(a*r_R*Z))*n;
        };
        return f;
    }
    std::function<double (const double *r)> getS_2(const Nucleus &nuc) const {
        double a = this->param;
        auto f = [a, nuc] (const double *r) -> double {
            const double *R = nuc.getCoord();
            double r_R = MathUtils::calcDistance(3, r, R);
            double Z = nuc.getCharge();
            double nomin = 2.0*(a - 1.0)*(1.0 - exp(a*r_R*Z)) - a*a*r_R*Z;
            double denom = 2.0*r_R*(1.0 + (a - 1.0)*exp(a*r_R*Z));
            return -Z*nomin/denom;
        };
        return f;
    }
protected:
    const double param;
};

#endif // SLATERCORRELATIONFUNCTION_H
