#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "RegularizedNuclearOperator.h"
#include "NuclearCorrelationFunction.h"
#include "Nucleus.h"
#include "utils/math_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

void RegularizedNuclearOperator::setupU_1(const NuclearCorrelationFunction &S, int d) {
    const Nuclei &nucs = this->nuclei;
    auto f = [nucs, &S, d] (const double *r) -> double {
        double result = 0.0;
        for (int i = 0; i < nucs.size(); i++) {
            result += S.getS_1(nucs[i], d)(r);
        }
        return result;
    };
    this->U_1[d].setReal(f);
}

void RegularizedNuclearOperator::setupU_2(const NuclearCorrelationFunction &S) {
    Printer::printHeader(0, "Setting up regularized nuclear potential");
    const Nuclei &nucs = this->nuclei;
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
    Printer::printSeparator(0, '=', 2);
}

} //namespace mrchem
