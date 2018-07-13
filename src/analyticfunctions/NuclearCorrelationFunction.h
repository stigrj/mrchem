#pragma once

#include "mrchem.h"
#include "utils/math_utils.h"

namespace mrchem {
class Nucleus;

class NuclearCorrelationFunction {
public:
    NuclearCorrelationFunction() = default;
    virtual ~NuclearCorrelationFunction() = default;

    virtual DoubleFunction getS_m1(const Nucleus &nuc) const = 0;
    virtual DoubleFunction getS_0(const Nucleus &nuc) const = 0;
    virtual DoubleFunction getS_1(const Nucleus &nuc, int d) const = 0;
    virtual DoubleFunction getS_2(const Nucleus &nuc) const = 0;
};

} //namespace mrchem
