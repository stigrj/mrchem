#pragma once

#include "NuclearOperator.h"
#include "NablaOperator.h"

namespace mrchem {
class NuclearCorrelationFunction;

class RegularizedNuclearOperator final : public NuclearOperator {
public:
    RegularizedNuclearOperator(const Nuclei &nucs,
                               const NuclearCorrelationFunction &S,
                               mrcpp::DerivativeOperator<3> &D)
            : NuclearOperator(nucs),
              nabla(D) {
        setupU_1(S, 0);
        setupU_1(S, 1);
        setupU_1(S, 2);
        setupU_2(S);

        RankZeroTensorOperator &h = (*this);
        h = -1.0*(U_1[0]*nabla[0] +
                  U_1[1]*nabla[1] +
                  U_1[2]*nabla[2] +
                  U_2);
    }

protected:
    NablaOperator nabla;
    AnalyticPotential U_1[3];

    void setupU_1(const NuclearCorrelationFunction &S, int d);
    void setupU_2(const NuclearCorrelationFunction &S);
};

} //namespace mrchem
