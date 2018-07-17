#pragma once

#include "QMPotential.h"
#include "RankZeroTensorOperator.h"

namespace mrchem {

class AnalyticPotential : public QMPotential {
public:
    AnalyticPotential(int adap = 1);
    virtual ~AnalyticPotential();

    void setReal(const DoubleFunction &func);
    void setImag(const DoubleFunction &func);

    virtual void setup(double prec);
    virtual void clear();

protected:
    mrcpp::AnalyticFunction<3> *real_func;
    mrcpp::AnalyticFunction<3> *imag_func;
};

class AnalyticPotentialOperator : public RankZeroTensorOperator {
public:
    AnalyticPotentialOperator() {
        RankZeroTensorOperator &v = (*this);
        v = pot;
    }

    void setReal(const DoubleFunction &func) { this->pot.setReal(func); }
    void setImag(const DoubleFunction &func) { this->pot.setImag(func); }


protected:
    AnalyticPotential pot;
};

} //namespace mrchem
