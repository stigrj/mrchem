#ifndef NUCLEARCORRELATIONOPERATOR_H
#define NUCLEARCORRELATIONOPERATOR_H

#include "AnalyticPotential.h"
#include "NuclearCorrelationFunction.h"
#include "Nucleus.h"

class NuclearCorrelationOperator : public AnalyticPotential {
public:
    NuclearCorrelationOperator(const Nuclei &nucs,
                               const NuclearCorrelationFunction &S) {
        const Nucleus &nuc = nucs[0];
        setReal(S.getS_0(nuc));
    }
    virtual ~NuclearCorrelationOperator() {
    }
};

#endif // NUCLEARCORRELATIONOPERATOR_H
