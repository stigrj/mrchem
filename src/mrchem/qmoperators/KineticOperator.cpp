#include "KineticOperator.h"
#include "OrbitalVector.h"
#include "Orbital.h"
#include "Timer.h"

using namespace std;
using namespace Eigen;

double KineticOperator::operator() (Orbital &orb_i, Orbital &orb_j, QMOperator *R) {
    RankZeroTensorOperator O_x;
    RankZeroTensorOperator O_y;
    RankZeroTensorOperator O_z;
    if (R != 0) {
        O_x = this->p[0]*(*R);
        O_y = this->p[1]*(*R);
        O_z = this->p[2]*(*R);
    } else {
        O_x = this->p[0];
        O_y = this->p[1];
        O_z = this->p[2];
    }

    Orbital *xOrb_j = O_x(orb_j);
    Orbital *xOrb_i = O_x(orb_i);
    complex<double> T_x = xOrb_i->dot(*xOrb_j);
    delete xOrb_i;
    delete xOrb_j;

    Orbital *yOrb_j = O_y(orb_j);
    Orbital *yOrb_i = O_y(orb_i);
    complex<double> T_y = yOrb_i->dot(*yOrb_j);
    delete yOrb_i;
    delete yOrb_j;

    Orbital *zOrb_j = O_z(orb_j);
    Orbital *zOrb_i = O_z(orb_i);
    complex<double> T_z = zOrb_i->dot(*zOrb_j);
    delete zOrb_i;
    delete zOrb_j;

    complex<double> T_tot = T_x + T_y + T_z;
    if (T_tot.imag() > MachineZero) {
        MSG_ERROR("Hermitian operator should have real expectation value");
    }
    return 0.5*T_tot.real();
}

double KineticOperator::adjoint(Orbital &orb_i, Orbital &orb_j, QMOperator *R) {
    NOT_IMPLEMENTED_ABORT;
}

MatrixXd KineticOperator::operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs, QMOperator *R) {
    Timer timer;
    TelePrompter::printHeader(1, "Compute Kinetic Matrix Elements");

    if (R != 0) return RankZeroTensorOperator::operator()(i_orbs, j_orbs, R);

    RankZeroTensorOperator &p_x = this->p[0];
    RankZeroTensorOperator &p_y = this->p[1];
    RankZeroTensorOperator &p_z = this->p[2];

    int Ni = i_orbs.size();
    int Nj = j_orbs.size();
    MatrixXcd T_x = MatrixXcd::Zero(Ni, Nj);
    MatrixXcd T_y = MatrixXcd::Zero(Ni, Nj);
    MatrixXcd T_z = MatrixXcd::Zero(Ni, Nj);
    {
        Timer timer;
        OrbitalVector dOrbs_j(0);
        for (int j = 0; j < Nj; j++) {
            Orbital &orb_j = j_orbs.getOrbital(j);
            Orbital *dOrb_j = p_x(orb_j);
            dOrbs_j.push_back(*dOrb_j);
        }
        for (int i = 0; i < Ni; i++) {
            Orbital &orb_i = i_orbs.getOrbital(i);
            Orbital *dOrb_i = p_x(orb_i);
            for (int j = 0; j < Nj; j++) {
                Orbital &dOrb_j = dOrbs_j.getOrbital(j);
                T_x(i,j) = dOrb_i->dot(dOrb_j);
            }
            delete dOrb_i;
        }
        dOrbs_j.clear(true);
        timer.stop();
        double t = timer.getWallTime();
        TelePrompter::printDouble(1, "T_x", t);
    }
    {
        Timer timer;
        OrbitalVector dOrbs_j(0);
        for (int j = 0; j < Nj; j++) {
            Orbital &orb_j = j_orbs.getOrbital(j);
            Orbital *dOrb_j = p_y(orb_j);
            dOrbs_j.push_back(*dOrb_j);
        }
        for (int i = 0; i < Ni; i++) {
            Orbital &orb_i = i_orbs.getOrbital(i);
            Orbital *dOrb_i = p_y(orb_i);
            for (int j = 0; j < Nj; j++) {
                Orbital &dOrb_j = dOrbs_j.getOrbital(j);
                T_y(i,j) = dOrb_i->dot(dOrb_j);
            }
            delete dOrb_i;
        }
        dOrbs_j.clear(true);
        timer.stop();
        double t = timer.getWallTime();
        TelePrompter::printDouble(1, "T_y", t);
    }
    {
        Timer timer;
        OrbitalVector dOrbs_j(0);
        for (int j = 0; j < Nj; j++) {
            Orbital &orb_j = j_orbs.getOrbital(j);
            Orbital *dOrb_j = p_z(orb_j);
            dOrbs_j.push_back(*dOrb_j);
        }
        for (int i = 0; i < Ni; i++) {
            Orbital &orb_i = i_orbs.getOrbital(i);
            Orbital *dOrb_i = p_z(orb_i);
            for (int j = 0; j < Nj; j++) {
                Orbital &dOrb_j = dOrbs_j.getOrbital(j);
                T_z(i,j) = dOrb_i->dot(dOrb_j);
            }
            delete dOrb_i;
        }
        dOrbs_j.clear(true);
        timer.stop();
        double t = timer.getWallTime();
        TelePrompter::printDouble(1, "T_z", t);
    }
    MatrixXcd T_tot = T_x + T_y + T_z;

    timer.stop();
    TelePrompter::printFooter(1, timer, 2);
    if (T_tot.imag().norm() > MachineZero) {
        MSG_ERROR("Hermitian operator should have real expectation value");
    }
    return 0.5*T_tot.real();
}

MatrixXd KineticOperator::adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs, QMOperator *R) {
    NOT_IMPLEMENTED_ABORT;
}

