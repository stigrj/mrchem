#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "AnalyticPotential.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

AnalyticPotential::AnalyticPotential(const DoubleFunction *f_r,
                                     const DoubleFunction *f_i)
        : QMPotential(1),
          real_func(nullptr),
          imag_func(nullptr) {
    if (f_r != nullptr) setReal(*f_r);
    if (f_i != nullptr) setImag(*f_i);
}

AnalyticPotential::~AnalyticPotential() {
    if (this->real_func != 0) delete this->real_func;
    if (this->imag_func != 0) delete this->imag_func;
}

void AnalyticPotential::setReal(const DoubleFunction &func) {
    if (this->real_func != 0) MSG_ERROR("Function not empty");
    this->real_func = new mrcpp::AnalyticFunction<3>(func);
}

void AnalyticPotential::setImag(const DoubleFunction &func) {
    if (this->imag_func != 0) MSG_ERROR("Function not empty");
    this->imag_func = new mrcpp::AnalyticFunction<3>(func);
}

void AnalyticPotential::setup(double prec) {
    if (this->isSetup(prec)) return;
    this->setApplyPrec(prec);

    if (this->hasReal()) MSG_ERROR("Potential not properly cleared");
    if (this->hasImag()) MSG_ERROR("Potential not properly cleared");

    Timer timer;
    if (this->real_func != 0) {
        alloc(NUMBER::Real);
        mrcpp::project(this->apply_prec, this->real(), *this->real_func);
    }
    if (this->imag_func != 0) {
        alloc(NUMBER::Imag);
        mrcpp::project(this->apply_prec, this->imag(), *this->imag_func);
    }
    timer.stop();

    int n = this->getNNodes();
    double t = timer.getWallTime();
    Printer::printTree(1, "Analytic potential", n, t);
}

void AnalyticPotential::clear() {
    this->free();           // delete FunctionTree pointers
    this->clearApplyPrec(); // apply_prec = -1
}

} //namespace mrchem
