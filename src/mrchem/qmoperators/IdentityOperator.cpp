#include "IdentityOperator.h"
#include "QMTensorOperator.h"
#include "OrbitalVector.h"
#include "Orbital.h"

using namespace Eigen;

extern MultiResolutionAnalysis<3> *MRA; // Global MRA
extern OrbitalVector workOrbVec;

IdentityOperator::IdentityOperator() : QMOperator(MRA->getMaxScale()) {
}

Orbital* IdentityOperator::operator() (Orbital &phi) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    Orbital *newOrb = new Orbital(phi);
    newOrb->deepCopy(phi);
    return newOrb;
}

Orbital* IdentityOperator::adjoint(Orbital &phi) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    NOT_IMPLEMENTED_ABORT;
}

double IdentityOperator::operator() (Orbital &phi_i, Orbital &phi_j, QMOperator *R) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    NOT_IMPLEMENTED_ABORT;
}

double IdentityOperator::adjoint(Orbital &phi_i, Orbital &phi_j, QMOperator *R) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    NOT_IMPLEMENTED_ABORT;
}

MatrixXd IdentityOperator::operator() (OrbitalVector &i_orbs, OrbitalVector &j_orbs, QMOperator *R) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    MatrixXcd S = MatrixXcd::Zero(i_orbs.size(), j_orbs.size());

    if (mpiOrbSize > 1) {
        if (&i_orbs == &j_orbs) {
            S = calcOverlapMatrix_P_H(i_orbs, j_orbs, R);
        } else {
            S = calcOverlapMatrix_P(i_orbs, j_orbs, R);
        }
    } else {
        S = calcOverlapMatrix(i_orbs, j_orbs, R);
    }

    if (S.imag().norm() > MachineZero) MSG_ERROR("Cannot handle complex yet");
    return S.real();
}

MatrixXd IdentityOperator::adjoint(OrbitalVector &i_orbs, OrbitalVector &j_orbs, QMOperator *R) {
    if (this->apply_prec < 0.0) MSG_ERROR("Uninitialized operator");

    NOT_IMPLEMENTED_ABORT;
}

MatrixXcd IdentityOperator::calcOverlapMatrix(OrbitalVector &bra, OrbitalVector &ket, QMOperator *R) {
    int Ni = bra.size();
    int Nj = ket.size();
    MatrixXcd S = MatrixXcd::Zero(Ni, Nj);

    if (R != 0) {
        RankZeroTensorOperator RR = (*R) * (*R);
        for (int j = 0; j < Nj; j++) {
            Orbital &phi_j = ket.getOrbital(j);
            Orbital *ket_j = RR(phi_j);
            for (int i = 0; i < Ni; i++) {
                Orbital &bra_i = bra.getOrbital(i);
                S(i,j) = bra_i.dot(*ket_j);
            }
            delete ket_j;
        }
    } else {
        for (int j = 0; j < Nj; j++) {
            Orbital &ket_j = ket.getOrbital(j);
            for (int i = 0; i < Ni; i++) {
                Orbital &bra_i = bra.getOrbital(i);
                S(i,j) = bra_i.dot(ket_j);
            }
        }
    }
    return S;
}

/** Calculate overlap matrix between two orbital sets */
MatrixXcd IdentityOperator::calcOverlapMatrix_P(OrbitalVector &bra, OrbitalVector &ket, QMOperator *R) {
#ifdef HAVE_MPI
    int Ni = bra.size();
    int Nj = ket.size();
    MatrixXcd S = MatrixXcd::Zero(Ni, Nj);

    OrbitalVector orbVecChunk_i(0);	//to store adresses of own i_orbs
    OrbitalVector orbVecChunk_j(0);	//to store adresses of own j_orbs
    OrbitalVector rcvOrbs(0);		//to store adresses of received orbitals

    vector<int> orbsIx;                 //to store own orbital indices
    int rcvOrbsIx[workOrbVecSize];	//to store received orbital indices

    //make vector with adresses of own orbitals
    for (int ix = mpiOrbRank; ix < Ni; ix += mpiOrbSize) {
        orbVecChunk_i.push_back(bra.getOrbital(ix));//i orbitals
        orbsIx.push_back(ix);
    }
    for (int jx = mpiOrbRank; jx < Nj; jx += mpiOrbSize) {
        orbVecChunk_j.push_back(ket.getOrbital(jx));//j orbitals
    }
    for (int iter = 0; iter >= 0; iter++) {
        //get a new chunk from other processes
        orbVecChunk_i.getOrbVecChunk(orbsIx, rcvOrbs, rcvOrbsIx, Ni, iter);

        if (R != 0) {
            //overlap between i and j chunks
            RankZeroTensorOperator RR = (*R) * (*R);
            for (int j = 0; j < orbVecChunk_j.size(); j++) {
                int jx = mpiOrbRank + j*mpiOrbSize;
                Orbital &phi_j = orbVecChunk_j.getOrbital(j);
                Orbital *ket_j = RR(phi_j);
                for (int i = 0; i < rcvOrbs.size(); i++) {
                    int ix = rcvOrbsIx[i];
                    Orbital &bra_i = rcvOrbs.getOrbital(i);
                    S(ix,jx) = bra_i.dot(*ket_j);
                }
                delete ket_i;
            }
        } else {
            //overlap between i and j chunks
            for (int j = 0; j < orbVecChunk_j.size(); j++) {
                int jx = mpiOrbRank + j*mpiOrbSize;
                Orbital &ket_j = orbVecChunk_j.getOrbital(j);
                for (int i = 0; i < rcvOrbs.size(); i++) {
                    int ix = rcvOrbsIx[i];
                    Orbital &bra_i = rcvOrbs.getOrbital(i);
                    S(ix,jx) = bra_i.dot(ket_j);
                }
            }
        }
        rcvOrbs.clearVec(false);//reset to zero size orbital vector
    }

    //clear orbital adresses (not the orbitals)
    orbVecChunk_i.clearVec(false);
    orbVecChunk_j.clearVec(false);
    workOrbVec.clear();

    MPI_Allreduce(MPI_IN_PLACE, &S(0,0), Ni*Nj,
                  MPI_DOUBLE_COMPLEX, MPI_SUM, mpiCommOrb);

    return S;
#else
    NOT_REACHED_ABORT;
#endif
}

/** Calculate overlap matrix between two orbital sets using MPI
 * 	assumes Hermitian overlap
 */
MatrixXcd IdentityOperator::calcOverlapMatrix_P_H(OrbitalVector &bra, OrbitalVector &ket, QMOperator *R) {
#ifdef HAVE_MPI
    int Ni = bra.size();
    int Nj = ket.size();
    MatrixXcd S = MatrixXcd::Zero(Ni, Nj);

    OrbitalVector orbVecChunk_i(0);	//to store adresses of own i_orbs
    OrbitalVector orbVecChunk_j(0);	//to store adresses of own j_orbs
    OrbitalVector rcvOrbs(0);		//to store adresses of received orbitals

    vector<int> orbsIx;                 //to store own orbital indices
    int rcvOrbsIx[workOrbVecSize];	//to store received orbital indices

    //make vector with adresses of own orbitals
    for (int ix = mpiOrbRank; ix < Ni; ix += mpiOrbSize) {
        orbVecChunk_i.push_back(bra.getOrbital(ix));//i orbitals
        orbsIx.push_back(ix);
    }
    for (int jx = mpiOrbRank; jx < Nj; jx += mpiOrbSize) {
        orbVecChunk_j.push_back(ket.getOrbital(jx));//j orbitals
    }

    //NB: last iteration may give empty chunk
    for (int iter = 0; iter >= 0; iter++) {
        //get a new chunk from other processes
        orbVecChunk_i.getOrbVecChunk_sym(orbsIx, rcvOrbs, rcvOrbsIx, Ni, iter);

        if (R != 0) {
            RankZeroTensorOperator RR = (*R) * (*R);
            //compute overlap between chunks
            for (int i = 0; i < rcvOrbs.size(); i++) {
                int ix = rcvOrbsIx[i];
                Orbital &phi_i = rcvOrbs.getOrbital(i);
                Orbital *bra_i = RR(phi_i);
                for (int j = 0; j < orbVecChunk_j.size(); j++) {
                    int jx = mpiOrbRank + j*mpiOrbSize;
                    if (ix%mpiOrbSize != mpiOrbRank or jx <= rcvOrbsIx[i]) {
                        Orbital &ket_j = orbVecChunk_j.getOrbital(j);
                        //compute only lower part in own block
                        S(ix,jx) = bra_i->dot(ket_j);
                        S(jx,ix) = conj(S(ix,jx));//symmetric
                    }
                }
                delete bra_i;
            }
        } else {
            //compute overlap between chunks
            for (int i = 0; i < rcvOrbs.size(); i++) {
                int ix = rcvOrbsIx[i];
                Orbital &bra_i = rcvOrbs.getOrbital(i);
                for (int j = 0; j < orbVecChunk_j.size(); j++) {
                    int jx = mpiOrbRank + j*mpiOrbSize;
                    //compute only lower part in own block
                    if (ix%mpiOrbSize != mpiOrbRank or jx <= rcvOrbsIx[i]) {
                        Orbital &ket_j = orbVecChunk_j.getOrbital(j);
                        S(ix,jx) = bra_i.dot(ket_j);
                        S(jx,ix) = conj(S(ix,jx));//symmetric
                    }
                }
            }
        }
        rcvOrbs.clearVec(false);//reset to zero size orbital vector
    }

    //clear orbital adresses (not the orbitals)
    orbVecChunk_i.clearVec(false);
    orbVecChunk_j.clearVec(false);
    workOrbVec.clear();

    MPI_Allreduce(MPI_IN_PLACE, &S(0,0), Ni*Nj,
                  MPI_DOUBLE_COMPLEX, MPI_SUM, mpiCommOrb);

    return S;
#else
    NOT_REACHED_ABORT;
#endif
}
