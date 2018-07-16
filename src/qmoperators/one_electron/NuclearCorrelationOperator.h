#pragma once

#include "RankZeroTensorOperator.h"
#include "AnalyticPotential.h"
#include "NuclearCorrelationFunction.h"
#include "Nucleus.h"

namespace mrchem {

class NuclearCorrelationOperator final : public RankZeroTensorOperator {
public:
    NuclearCorrelationOperator(const Nuclei &nucs, const NuclearCorrelationFunction &S) {
        std::vector<DoubleFunction> funcs;
        for (int k = 0; k < nucs.size(); k++) {
            funcs.push_back(S.getS_0(nucs[k]));
        }

        auto f1 = [funcs] (const double *r) -> double {
            double result = 1.0;
            for (int i = 0; i < funcs.size(); i++) {
                result *= funcs[i](r);
            }
            return result;
        };
        this->R.setReal(f1);

        auto f2 = [funcs] (const double *r) -> double {
            double result = 1.0;
            for (int i = 0; i < funcs.size(); i++) {
                result *= std::pow(funcs[i](r), 2.0);
            }
            return result;
        };
        this->R2.setReal(f2);

        RankZeroTensorOperator &V = (*this);
        V = R;
    }

    RankZeroTensorOperator getR2() { return this->R*this->R; }

    void setup(double prec) {
        RankZeroTensorOperator::setup(prec);
        this->R2.setup(prec);
    }

    void clear() {
        RankZeroTensorOperator::clear();
        this->R2.clear();
    }

protected:
    AnalyticPotential R;
    AnalyticPotential R2;
};

} //namespace mrchem
