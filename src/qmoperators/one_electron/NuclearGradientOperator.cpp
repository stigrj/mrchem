#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "NuclearGradientOperator.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

void NuclearGradientPotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);

    if (hasReal()) MSG_ERROR("Potential not properly cleared");
    if (hasImag()) MSG_ERROR("Potential not properly cleared");

    Timer timer;
    alloc(NUMBER::Real);
    mrcpp::build_grid(this->real(), this->func);
    mrcpp::project(this->apply_prec, this->real(), this->func);
    timer.stop();

    int n = getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(0, "Nuclear potential", n, t);
}

void NuclearGradientPotential::clear() {
    free();           // delete FunctionTree pointers
    clearApplyPrec(); // apply_prec = -1
}

} //namespace mrchem
