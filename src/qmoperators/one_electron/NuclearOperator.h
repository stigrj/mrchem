#pragma once

#include "RankZeroTensorOperator.h"
#include "AnalyticPotential.h"
#include "Nucleus.h"

namespace mrchem {

class NuclearOperator : public RankZeroTensorOperator {
public:
    NuclearOperator(const Nuclei &nucs) : nuclei(nucs) { }
    virtual ~NuclearOperator() = 0; // To avoid instantiation of base class

    Nuclei &getNuclei() { return this->nuclei; }
    const Nuclei &getNuclei() const { return this->nuclei; }

    double trace(const Nuclei &nucs_b) const;
    using RankZeroTensorOperator::trace;

protected:
    Nuclei nuclei;
    AnalyticPotential U_2;
};

} //namespace mrchem
