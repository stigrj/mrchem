#pragma once

#pragma GCC system_header
#include <Eigen/Core>

#include "QMOperator.h"

class IdentityOperator : public QMOperator {
public:
    IdentityOperator();
    virtual ~IdentityOperator() { }

    virtual void setup(double prec) { setApplyPrec(prec); }
    virtual void clear() { clearApplyPrec(); }

    virtual Orbital* operator() (Orbital &phi);
    virtual Orbital* adjoint(Orbital &phi);

    virtual double operator() (Orbital &phi_i, Orbital &phi_j, QMOperator *R = 0);
    virtual double adjoint(Orbital &phi_i, Orbital &phi_j, QMOperator *R = 0);

    virtual Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs, QMOperator *R = 0);
    virtual Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs, QMOperator *R = 0);

protected:
    Eigen::MatrixXcd calcOverlapMatrix(OrbitalVector &bra, OrbitalVector &ket, QMOperator *R = 0);
    Eigen::MatrixXcd calcOverlapMatrix_P(OrbitalVector &bra, OrbitalVector &ket, QMOperator *R = 0);
    Eigen::MatrixXcd calcOverlapMatrix_P_H(OrbitalVector &bra, OrbitalVector &ket, QMOperator *R = 0);
};

