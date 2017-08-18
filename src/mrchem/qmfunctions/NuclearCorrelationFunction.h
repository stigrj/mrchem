#ifndef NUCLEARCORRELATIONFUNCTION_H
#define NUCLEARCORRELATIONFUNCTION_H

class Nucleus;

class NuclearCorrelationFunction {
public:
    NuclearCorrelationFunction() { }
    virtual ~NuclearCorrelationFunction() { }

    virtual std::function<double (const double *r)> getS_m1(const Nucleus &nuc) const = 0;
    virtual std::function<double (const double *r)> getS_0(const Nucleus &nuc) const = 0;
    virtual std::function<double (const double *r)> getS_1(const Nucleus &nuc, int d) const = 0;
    virtual std::function<double (const double *r)> getS_2(const Nucleus &nuc) const = 0;
};

#endif // NUCLEARCORRELATIONFUNCTION_H
