#include "CoulombPotential.h"
#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "parallel.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/orbital_utils.h"

using mrcpp::PoissonOperator;
using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief constructor
 *
 * @param P: MW Poisson operator
 * @param Phi: orbitals defining the operator
 *
 * Density is not spin-separated since this operator requires only total density.
 * This operator will always point to the same OrbitalVector, but the orbitals within
 * the vector can change throughout the calculation. The density and (*this)
 * QMPotential is uninitialized at this point and will be computed at setup.
 */

CoulombPotential::CoulombPotential(PoissonOperator *P, OrbitalVector *Phi)
        : QMPotential(1, mpi::share_coul_pot)
        , local(mpi::local_coul_pot)
        , density(mpi::share_coul_dens)
        , orbitals(Phi)
        , poisson(P) {}

/** @brief prepare operator for application
 *
 * @param prec: apply precision
 *
 * This will compute the Coulomb potential by application of the Poisson
 * operator to the density. If the density is not available it is computed
 * from the current orbitals (assuming that the orbitals are available).
 * For first-order perturbations the first order density and the Hessian will be
 * computed. In order to make the Hessian available to CoulombOperator, it is stored in the 
 * potential function instead of the zeroth-order potential.
 *
 */
void CoulombPotential::setup(double prec) {
    if (isSetup(prec)) return;
    setApplyPrec(prec);
    if (useLocal()) {
        setupLocalDensity(prec);
        QMFunction V = setupLocalPotential(prec);
        allreducePotential(prec, V);
    } else {
        setupDensity(prec);
        setupPotential(prec);
    }
}

/** @brief clear operator after application
 *
 * This will clear the operator and bring it back to the state after construction.
 * The operator can now be reused after another setup.
 */
void CoulombPotential::clear() {
    QMFunction::free(NUMBER::Total);   // delete FunctionTree pointers
    this->density.free(NUMBER::Total); // delete FunctionTree pointers
    clearApplyPrec();                  // apply_prec = -1
}

/** @brief compute Coulomb potential
 *
 * @param prec: apply precision
 *
 * This will compute the Coulomb potential by application o the Poisson operator
 * to the precomputed electron density.
 */
void CoulombPotential::setupPotential(double prec) {
    if (this->poisson == nullptr) MSG_ERROR("Poisson operator not initialized");

    PoissonOperator &P = *this->poisson;
    QMFunction &V = *this;
    QMFunction &rho = this->density;

    if (V.hasReal()) MSG_ERROR("Potential not properly cleared");
    if (V.hasImag()) MSG_ERROR("Potential not properly cleared");

    // Adjust precision by system size
    double abs_prec = prec / rho.norm();
    bool need_to_apply = not(V.isShared()) or mpi::share_master();

    Timer timer;
    V.alloc(NUMBER::Real);
    if (need_to_apply) mrcpp::apply(abs_prec, V.real(), P, rho.real());
    mpi::share_function(V, 0, 22445, mpi::comm_share);
    timer.stop();

    int n = V.getNNodes(NUMBER::Total);
    double t = timer.getWallTime();
    Printer::printTree(0, "Coulomb potential", n, t);
}

/** @brief compute Coulomb potential
 *
 * @param prec: apply precision
 *
 * This will compute the Coulomb potential by application o the Poisson operator
 * to the precomputed electron density.
 */
QMFunction CoulombPotential::setupLocalPotential(double prec) {
    if (this->poisson == nullptr) MSG_ERROR("Poisson operator not initialized");

    PoissonOperator &P = *this->poisson;
    OrbitalVector &Phi = *this->orbitals;
    QMFunction &rho = this->density;

    // Adjust precision by system size
    double abs_prec = prec / orbital::get_electron_number(Phi);

    Timer timer;
    QMFunction V(false);
    V.alloc(NUMBER::Real);
    mrcpp::apply(abs_prec, V.real(), P, rho.real());
    timer.stop();

    int n = V.getNNodes(NUMBER::Total);
    double t = timer.getWallTime();
    Printer::printTree(0, "Coulomb potential", n, t);

    return V;
}

void CoulombPotential::allreducePotential(double prec, QMFunction &V_loc) {
    Timer t_com;

    QMFunction &V_tot = *this;
    OrbitalVector &Phi = *this->orbitals;

    double abs_prec = prec / orbital::get_electron_number(Phi);

    // Add up local contributions into the grand master
    mpi::reduce_function(abs_prec, V_loc, mpi::comm_orb);

    if (not V_tot.hasReal()) V_tot.alloc(NUMBER::Real);
    if (V_tot.isShared()) {
        int tag = 3141;
        // MPI grand master distributes to shared masters
        mpi::broadcast_function(V_loc, mpi::comm_sh_group);
        if (mpi::share_master()) {
            // MPI shared masters copies the function into final memory
            mrcpp::copy_grid(V_tot.real(), V_loc.real());
            mrcpp::copy_func(V_tot.real(), V_loc.real());
        }
        // MPI share masters distributes to their sharing ranks
        mpi::share_function(V_tot, 0, tag, mpi::comm_share);
    } else {
        // MPI grand master distributes to all ranks
        mpi::broadcast_function(V_loc, mpi::comm_orb);
        // All MPI ranks copies the function into final memory
        mrcpp::copy_grid(V_tot.real(), V_loc.real());
        mrcpp::copy_func(V_tot.real(), V_loc.real());
    }
    t_com.stop();
    Printer::printTree(0, "Allreduce potential", V_tot.getNNodes(NUMBER::Total), t_com.getWallTime());
}

} //namespace mrchem
