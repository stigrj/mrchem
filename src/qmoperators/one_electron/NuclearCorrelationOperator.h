#pragma once

#include "QMPotential.h"
#include "RankZeroTensorOperator.h"

namespace mrchem {
class NuclearCorrelationFunction;

class NuclearCorrelationPotential final : public QMPotential {
public:
    NuclearCorrelationPotential(const Nuclei &nucs, const NuclearCorrelationFunction &S);

protected:
    mrcpp::AnalyticFunction<3> func;

    void setup(double prec);
    void clear();
};

class NuclearCorrelationOperator final : public RankZeroTensorOperator {
public:
    NuclearCorrelationOperator(const Nuclei &nucs, const NuclearCorrelationFunction &S)
            : ncp(nucs, S) {
        RankZeroTensorOperator &R = (*this);
        R = ncp;
    }

protected:
    NuclearCorrelationPotential ncp;
};

} //namespace mrchem
