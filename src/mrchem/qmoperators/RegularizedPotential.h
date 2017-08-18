#ifndef REGULARIZEDPOTENTIAL_H
#define REGULARIZEDPOTENTIAL_H

#include "QMTensorOperator.h"
#include "AnalyticPotential.h"
#include "NablaOperator.h"
#include "NuclearCorrelationFunction.h"
#include "Nucleus.h"

class RegularizedPotential : public RankZeroTensorOperator {
public:
    RegularizedPotential(Nuclei &nucs,
                         NuclearCorrelationFunction &R,
                         DerivativeOperator<3> &D)
            : nuclei(nucs), nabla(D) {
        const Nucleus &nuc = this->nuclei[0];
        U_x.setReal(R.getS_1(nuc, 0));
        U_y.setReal(R.getS_1(nuc, 1));
        U_z.setReal(R.getS_1(nuc, 2));
        U_2.setReal(R.getS_2(nuc));
        initializeTensorOperator();
    }
    virtual ~RegularizedPotential() { }

    Nuclei &getNuclei() { return this->nuclei; }
    const Nuclei &getNuclei() const { return this->nuclei; }

protected:
    Nuclei nuclei;
    NablaOperator nabla;
    AnalyticPotential U_x;
    AnalyticPotential U_y;
    AnalyticPotential U_z;
    AnalyticPotential U_2;

    void initializeTensorOperator() {
        RankZeroTensorOperator &h = (*this);
        h = -1.0*(U_x*nabla[0] +
                  U_y*nabla[1] +
                  U_z*nabla[2] +
                  U_2);
    }
};

#endif // REGULARIZEDPOTENTIAL_H
