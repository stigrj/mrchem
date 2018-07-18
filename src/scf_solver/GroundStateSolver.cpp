#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"

#include "GroundStateSolver.h"
#include "FockOperator.h"
#include "Orbital.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief constructor
 *
 * @param h: Helmholtz operators
 *
 * SCF solver will NOT take ownership of the HelmholtzVector, so the original object
 * must be taken care of externally (do not delete until SCF goes out of scope).
 * Fock matrix, FockOperator and OrbitalVector are not initialized at this stage,
 * so the SCF solver needs to be "setup()" before "optimize()".
 */
GroundStateSolver::GroundStateSolver(HelmholtzVector &h)
        : SCF(h),
          fMat_n(nullptr),
          fOper_n(nullptr),
          orbitals_n(nullptr),
          nuc_corr_fac(nullptr) {
}

/** @brief Computes the Helmholtz argument for the all orbitals.
 *
 * @param fock: Fock operator (potential part used in first term)
 * @param M: Rotation matrix for second term
 * @param Phi: Orbital vector
 * @param clearFock: clear Fock operator after application
 *
 * Argument contains the potential operator acting on orbital i, and the sum
 * of all orbitals weighted by the Fock matrix. The effect of using inexact
 * Helmholtz operators are included in Lambda, which is a diagonal matrix
 * with the actual lambda parameters used in the Helmholtz operators
 * (input matrix M is assumed to be L-F).
 *
 * psi_j = \hat{V}phi_i + \sum_j (\Lambda_{ij}-F_{ij})phi_j
 *
 */
OrbitalVector GroundStateSolver::setupHelmholtzArguments(FockOperator &fock,
                                                         const ComplexMatrix &M,
                                                         OrbitalVector &Phi,
                                                         bool clearFock) {
    Timer timer_tot;
    Printer::printHeader(0, "Setting up Helmholtz arguments");
    int oldprec = Printer::setPrecision(5);

    RankZeroTensorOperator &V = fock.potential();
    
    double coef = -1.0/(2.0*MATHCONST::pi);
    Timer timer_1;
    OrbitalVector part_1 = orbital::param_copy(Phi);
    for (int i = 0; i < Phi.size(); i++) {
        if (mpi::my_orb(Phi[i])) {
            part_1[i] = V(Phi[i]);
        }
    }
    timer_1.stop();
    Printer::printDouble(0, "Potential part", timer_1.getWallTime(), 5);

    if (clearFock) fock.clear();

    Timer timer_2;
    OrbitalVector part_2 = orbital::multiply(M, Phi);
    timer_2.stop();
    Printer::printDouble(0, "Matrix part", timer_2.getWallTime(), 5);

    OrbitalVector out = orbital::add(coef, part_1, coef, part_2);
    orbital::free(part_1);
    orbital::free(part_2);

    /*
    Timer timer_2;
    OrbitalVector orbVecChunk_i(0); //to store adresses of own i_orbs
    OrbitalVector rcvOrbs(0);       //to store adresses of received orbitals
    vector<int> orbsIx;             //to store own orbital indices
    int rcvOrbsIx[workOrbVecSize];  //to store received orbital indices

    //make vector with adresses of own orbitals
    int Ni = phi.size();
    int maxOrbPerMpi = Ni/mpiOrbSize + 1;//upper bound
    for (int ix = mpiOrbRank; ix < Ni; ix += mpiOrbSize) {
        orbVecChunk_i.push_back(phi.getOrbital(ix));//i orbitals
        orbsIx.push_back(ix);
    }

    for (int iter = 0; iter >= 0; iter++) {
        //get a new chunk from other processes
        orbVecChunk_i.getOrbVecChunk(orbsIx, rcvOrbs, rcvOrbsIx, Ni, iter);
        for (int i = mpiOrbRank; i < Ni; i += mpiOrbSize) {
            Orbital &phi_i = phi.getOrbital(i);

            vector<complex<double> > coefs;
            vector<Orbital *> orbs;

            int ix = i;
            for (int j = 0; j < rcvOrbs.size(); j++) {
                int jx = rcvOrbsIx[j];
                double coef = M(ix,jx);
                // Linear scaling screening inserted here
                // if (std::abs(coef) > MachineZero) {
                    Orbital &phi_j = rcvOrbs.getOrbital(j);
                    double norm_j = sqrt(phi_j.getSquareNorm());
                    if (norm_j > 0.01*getOrbitalPrecision()) {
                        coefs.push_back(coef);
                        orbs.push_back(&phi_j);
                    }
                //}
            }

            Orbital *tmp_i = new Orbital(phi_i);
            if (orbs.size() > 0) this->add(*tmp_i, coefs, orbs, false);

            this->add.inPlace(out[i], coef, *tmp_i);
            delete tmp_i;
        }
        rcvOrbs.clearVec(false);//reset to zero size orbital vector
    }
    orbVecChunk_i.clearVec(false);
    workOrbVec.clear();
    timer_2.stop();
    Printer::printDouble(0, "Matrix part", timer_2.getWallTime(), 5);
    */

    timer_tot.stop();
    Printer::printFooter(0, timer_tot, 2);
    Printer::setPrecision(oldprec);

    return out;
}

/** @brief Computes the SCF energy by tracing the Fock operator
 *
 * Prints the current nuclear, electronic and total energies.
 */
double GroundStateSolver::calcProperty() {
    Printer::printHeader(0, "Calculating SCF energy");
    Timer timer;

    ComplexMatrix &F = *this->fMat_n;
    FockOperator &fock = *this->fOper_n;
    OrbitalVector &Phi = *this->orbitals_n;
    NuclearCorrelationOperator *R = this->nuc_corr_fac;

    SCFEnergy E = fock.trace(Phi, F, R);
    this->energy.push_back(E);

    timer.stop();
    int oldPrec = Printer::setPrecision(15);
    println(0, " Nuclear energy              " << std::setw(30) << E.getNuclearEnergy());
    println(0, " Electronic energy           " << std::setw(30) << E.getElectronicEnergy());
    Printer::printSeparator(0, '-');
    println(0, " Total energy                " << std::setw(30) << E.getTotalEnergy());
    Printer::printFooter(0, timer, 2);
    Printer::setPrecision(oldPrec);

    return E.getTotalEnergy();
}

/** @brief Computes the SCF energy update from last iteration */
double GroundStateSolver::calcPropertyError() const {
    int iter = this->property.size();
    return std::abs(getUpdate(this->property, iter, false));
}

/** @brief Pretty printing of the different contributions to the SCF energy */
void GroundStateSolver::printProperty() const {
    SCFEnergy scf_0, scf_1;
    int iter = this->energy.size();
    if (iter > 1) scf_0 = this->energy[iter - 2];
    if (iter > 0) scf_1 = this->energy[iter - 1];

    double phi_0 = scf_0.getOrbitalEnergy();
    double phi_1 = scf_1.getOrbitalEnergy();
    double T_0 = scf_0.getKineticEnergy();
    double T_1 = scf_1.getKineticEnergy();
    double V_0 = scf_0.getElectronNuclearEnergy();
    double V_1 = scf_1.getElectronNuclearEnergy();
    double J_0 = scf_0.getElectronElectronEnergy();
    double J_1 = scf_1.getElectronElectronEnergy();
    double K_0 = scf_0.getExchangeEnergy();
    double K_1 = scf_1.getExchangeEnergy();
    double XC_0 = scf_0.getExchangeCorrelationEnergy();
    double XC_1 = scf_1.getExchangeCorrelationEnergy();
    double E_0 = scf_0.getElectronicEnergy();
    double E_1 = scf_1.getElectronicEnergy();
    double N_0 = scf_0.getNuclearEnergy();
    double N_1 = scf_1.getNuclearEnergy();

    Printer::printHeader(0, "                    Energy                 Update      Done ");
    printUpdate(" Orbital    ", phi_1, phi_1 - phi_0);
    printUpdate(" Kinetic    ",   T_1,   T_1 -   T_0);
    printUpdate(" N-E        ",   V_1,   V_1 -   V_0);
    printUpdate(" Coulomb    ",   J_1,   J_1 -   J_0);
    printUpdate(" Exchange   ",   K_1,   K_1 -   K_0);
    printUpdate(" X-C        ",  XC_1,  XC_1 -  XC_0);
    Printer::printSeparator(0, '-');
    printUpdate(" Electronic ",  E_1,  E_1 -  E_0);
    printUpdate(" Nuclear    ",  N_1,  N_1 -  N_0);
    Printer::printSeparator(0, '-');
    printUpdate(" Total      ",  E_1 + N_1, (E_1+N_1) - (E_0+N_0));
    Printer::printSeparator(0, '=');
}

} //namespace mrchem
