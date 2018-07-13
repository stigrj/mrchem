#pragma once

#include "RankZeroTensorOperator.h"
#include "AnalyticPotential.h"
#include "NuclearCorrelationFunction.h"

namespace mrchem {

class NuclearCorrelationOperator final : public RankZeroTensorOperator {
public:
    NuclearCorrelationOperator(const Nuclei &nucs, const NuclearCorrelationFunction &S) {
        std::vector<DoubleFunction> funcs;
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
        this->ncp.setReal(f);

        RankZeroTensorOperator &R = (*this);
        R = ncp;
    }

protected:
    AnalyticPotential ncp;
};

} //namespace mrchem
