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
                         NuclearCorrelationFunction &S,
                         DerivativeOperator<3> &D)
            : nuclei(nucs), nabla(D) {
        setupU_1(nucs, S, 0);
        setupU_1(nucs, S, 1);
        setupU_1(nucs, S, 2);
        setupU_2(nucs, S);
        initializeTensorOperator();
    }
    virtual ~RegularizedPotential() { }

    Nuclei &getNuclei() { return this->nuclei; }
    const Nuclei &getNuclei() const { return this->nuclei; }

protected:
    Nuclei nuclei;
    NablaOperator nabla;
    AnalyticPotential U_1[3];
    AnalyticPotential U_2;

    void initializeTensorOperator() {
        RankZeroTensorOperator &h = (*this);
        h = -1.0*(U_1[0]*nabla[0] +
                  U_1[1]*nabla[1] +
                  U_1[2]*nabla[2] +
                  U_2);
    }

    void setupU_1(const Nuclei &nucs, NuclearCorrelationFunction &S, int d) {
        auto f = [nucs, &S, d] (const double *r) -> double {
            double result = 0.0;
            for (int i = 0; i < nucs.size(); i++) {
                result += S.getS_1(nucs[i], d)(r);
            }
            return result;
        };

        this->U_1[d].setReal(f);
    }

    void setupU_2(const Nuclei &nucs, NuclearCorrelationFunction &S) {
        auto f = [nucs, &S] (const double *r) -> double {
            double result = 0.0;
            for (int a = 0; a < nucs.size(); a++) {
                const Nucleus &A = nucs[a];
                for (int b = 0; b < a; b++) {
                    const Nucleus &B = nucs[b];
                    result += S.getS_1(A, 0)(r)*S.getS_1(B, 0)(r);
                    result += S.getS_1(A, 1)(r)*S.getS_1(B, 1)(r);
                    result += S.getS_1(A, 2)(r)*S.getS_1(B, 2)(r);
                }
                result += S.getS_2(A)(r);
            }
            return result;
        };

        this->U_2.setReal(f);
    }
};

#endif // REGULARIZEDPOTENTIAL_H
