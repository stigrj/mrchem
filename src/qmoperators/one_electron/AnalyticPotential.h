#pragma once

#include "QMPotential.h"

namespace mrchem {

class AnalyticPotential : public QMPotential {
public:
    AnalyticPotential(const DoubleFunction *f_r = nullptr, const DoubleFunction *f_i = nullptr);
    virtual ~AnalyticPotential();

    void setReal(const DoubleFunction &func);
    void setImag(const DoubleFunction &func);

    virtual void setup(double prec);
    virtual void clear();

protected:
    mrcpp::AnalyticFunction<3> *real_func;
    mrcpp::AnalyticFunction<3> *imag_func;
};

} //namespace mrchem
