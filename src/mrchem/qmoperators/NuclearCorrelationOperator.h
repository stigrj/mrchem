#ifndef NUCLEARCORRELATIONOPERATOR_H
#define NUCLEARCORRELATIONOPERATOR_H

#include <vector>

#include "AnalyticPotential.h"
#include "NuclearCorrelationFunction.h"
#include "Nucleus.h"

class NuclearCorrelationOperator : public AnalyticPotential {
public:
    NuclearCorrelationOperator(const Nuclei &nucs,
                               const NuclearCorrelationFunction &S) {
        std::vector<std::function<double (const double *r)> > funcs;
        for (int k = 0; k < nucs.size(); k++) {
            funcs.push_back(S.getS_0(nucs[k]));
        }

        auto f = [funcs] (const double *r) -> double {
            double result = 1.0;
            for (int i = 0; i < funcs.size(); i++) {
                result *= funcs[i](r);
            }
            return result;
        };
        setReal(f);
    }
    virtual ~NuclearCorrelationOperator() {
    }
};

#endif // NUCLEARCORRELATIONOPERATOR_H
