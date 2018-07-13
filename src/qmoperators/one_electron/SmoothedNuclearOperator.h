#pragma once

#include "NuclearOperator.h"

namespace mrchem {

class SmoothedNuclearOperator final : public NuclearOperator {
public:
    SmoothedNuclearOperator(const Nuclei &nucs, double prec)
            : NuclearOperator(nucs) {
        setupU_2(prec);

        RankZeroTensorOperator &v = (*this);
        v = U_2;
    }

protected:
    void setupU_2(double prec);
};

} //namespace mrchem
