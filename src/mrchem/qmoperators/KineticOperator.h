#pragma once

#include "QMTensorOperator.h"
#include "MomentumOperator.h"

class KineticOperator : public RankZeroTensorOperator {
public:
    KineticOperator(DerivativeOperator<3> &D)
            : p(D) {
        initializeTensorOperator();
    }
    virtual ~KineticOperator() { }

    virtual double operator() (Orbital &orb_i, Orbital &orb_j, QMOperator *R = 0);
    virtual double adjoint(Orbital &orb_i, Orbital &orb_j, QMOperator *R = 0);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs, QMOperator *R = 0);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs, QMOperator *R = 0);

    using RankZeroTensorOperator::operator();
    using RankZeroTensorOperator::adjoint;

protected:
    MomentumOperator p;

    void initializeTensorOperator() {
        RankZeroTensorOperator &p_x = this->p[0];
        RankZeroTensorOperator &p_y = this->p[1];
        RankZeroTensorOperator &p_z = this->p[2];

        RankZeroTensorOperator &h = *this;
        h = -0.5*(p_x*p_x + p_y*p_y + p_z*p_z);
    }
};

