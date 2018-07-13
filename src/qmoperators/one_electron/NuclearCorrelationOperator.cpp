#include <vector>

#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "NuclearCorrelationOperator.h"
#include "NuclearCorrelationFunction.h"
#include "Nucleus.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

NuclearCorrelationPotential::NuclearCorrelationPotential(const Nuclei &nucs, const NuclearCorrelationFunction &S)
        : QMPotential(1) {
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
    this->func.set(f);
}

void NuclearCorrelationPotential::setup(double prec) {
    if (this->isSetup(prec)) return;
    this->setApplyPrec(prec);

    if (this->hasReal()) MSG_ERROR("Potential not properly cleared");
    if (this->hasImag()) MSG_ERROR("Potential not properly cleared");

    Timer timer;
    this->alloc(NUMBER::Real);
    mrcpp::project(this->apply_prec, this->real(), this->func);
    timer.stop();

    int n = this->getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(1, "NuclearCorrelationPotential", n, t);
}

void NuclearCorrelationPotential::clear() {
    this->free();           // delete FunctionTree pointers
    this->clearApplyPrec(); // apply_prec = -1
}

} //namespace mrchem
