#pragma once

#pragma GCC system_header
#include <Eigen/Core>
#include <vector>

#include "QMOperator.h"

class QMOperatorExp;
class RegularizedPotential;
class KineticOperator;
class CoulombOperator;
class ExchangeOperator;
class XCOperator;
class Orbital;
class OrbitalVector;
class SCFEnergy;

class FockOperator : public QMOperator {
public:
    FockOperator(KineticOperator *t = 0,
                 RegularizedPotential *v = 0,
                 CoulombOperator *j = 0,
                 ExchangeOperator *k = 0,
                 XCOperator *xc = 0);
    virtual ~FockOperator();

    void setPerturbationOperator(QMOperatorExp &h_1) { this->H_1 = &h_1; }
    QMOperatorExp *getPerturbationOperator() { return this->H_1; }

    KineticOperator *getKineticOperator() { return this->T; }
    RegularizedPotential *getNuclearPotential() { return this->V; }
    CoulombOperator *getCoulombOperator() { return this->J; }
    ExchangeOperator *getExchangeOperator() { return this->K; }
    XCOperator *getXCOperator() { return this->XC; }

    void rotate(Eigen::MatrixXd &U);

    virtual void setup(double prec);
    virtual void clear();

    Orbital* operator() (Orbital &orb_p);
    Orbital* adjoint(Orbital &orb_p);

    double operator() (Orbital &orb_i, Orbital &orb_j, QMOperator *R = 0);
    double adjoint(Orbital &orb_i, Orbital &orb_j, QMOperator *R = 0);

    Eigen::MatrixXd operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs, QMOperator *R = 0);
    Eigen::MatrixXd adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs, QMOperator *R = 0);

    Orbital* applyPotential(Orbital &orb_p);
    double applyPotential(Orbital &orb_i, Orbital &orb_j, QMOperator *R = 0);
    Eigen::MatrixXd applyPotential(OrbitalVector &i_orbs, OrbitalVector &j_orbs, QMOperator *R = 0);

    Orbital* applyAdjointPotential(Orbital &orb_p);
    double applyAdjointPotential(Orbital &orb_i, Orbital &orb_j, QMOperator *R = 0);
    Eigen::MatrixXd applyAdjointPotential(OrbitalVector &i_orbs, OrbitalVector &j_orbs, QMOperator *R = 0);

    SCFEnergy trace(OrbitalVector &phi, Eigen::MatrixXd &F, QMOperator *R = 0);

protected:
    KineticOperator *T;
    RegularizedPotential *V;
    CoulombOperator *J;
    ExchangeOperator *K;
    XCOperator *XC;
    QMOperatorExp *H_1;   // First order perturbation operators
};


