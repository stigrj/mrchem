#pragma once

#include "QMPotential.h"
#include "RankOneTensorOperator.h"
#include "NuclearGradientFunction.h"

namespace mrchem {

class NuclearGradientPotential final : public QMPotential {
public:
    NuclearGradientPotential(int d, const Nucleus &nuc, double c)
        : QMPotential(1), func(d, nuc, c) { }

    void setup(double prec);
    void clear();

protected:
    NuclearGradientFunction func;
};

class NuclearGradientOperator : public RankOneTensorOperator<3> {
public:
    NuclearGradientOperator(const Nucleus &nuc, double c)
            : x_rm3(0, nuc, c),
              y_rm3(1, nuc, c),
              z_rm3(2, nuc, c) {
        RankOneTensorOperator &v = (*this);
        v[0] = x_rm3;
        v[1] = y_rm3;
        v[2] = z_rm3;
    }

protected:
    NuclearGradientPotential x_rm3;
    NuclearGradientPotential y_rm3;
    NuclearGradientPotential z_rm3;
};


} //namespace mrchem
