#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include "utils/math_utils.h"

#include "NuclearOperator.h"
#include "parallel.h"
#include "qmfunctions/qmfunction_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief projects an analytic expression for a smoothed nuclear potential
 *
 * @param[in] nucs collection of nuclei that defines the potential
 * @param[in] prec precision used both in smoothing and projection
 *
 * We create two different analytic functions for the nuclear potential:
 * 
 * 1) this->func: The total potential from all nuclei of the system.
 *                This is needed later for analytic calculations.
 *
 * 2) loc_func: Temporary function that contains only some of the
 *              nuclei, which are distributed among the available
 *              MPIs. This is used only for the projection below.
 */
NuclearPotential::NuclearPotential(const Nuclei &nucs, double prec)
        : QMPotential(1, mpi::share_nuc_pot) {
    int oldprec = Printer::setPrecision(5);
    Printer::printHeader(0, "Setting up nuclear potential");
    println(0, " Nr  Element         Charge        Precision     Smoothing ");
    Printer::printSeparator(0, '-');

    NuclearFunction loc_func;

    double c = 0.00435 * prec;
    for (int i = 0; i < nucs.size(); i++) {
        const Nucleus &nuc = nucs[i];
        double Z = nuc.getCharge();
        double Z_5 = std::pow(Z, 5.0);
        double smooth = std::pow(c / Z_5, 1.0 / 3.0);

        this->func.push_back(nuc, smooth);
        if (mpi::orb_rank == i%mpi::orb_size) loc_func.push_back(nuc, smooth); 

        std::stringstream symbol;
        symbol << nuc.getElement().getSymbol();
        symbol << "  ";
        printout(0, std::setw(3) << i + 1 << "     ");
        printout(0, symbol.str()[0] << symbol.str()[1]);
        printout(0, std::setw(21) << Z);
        printout(0, std::setw(14) << prec);
        printout(0, std::setw(14) << smooth << std::endl);
    }

    Timer t_tot;

    Timer t_loc;
    QMFunction V_loc(false);
    qmfunction::project(V_loc, loc_func, NUMBER::Real, prec);
    t_loc.stop();
    Printer::printSeparator(0, '-');
    Printer::printTree(0, "Nuclear potential", V_loc.getNNodes(NUMBER::Total), t_loc.getWallTime());

    Timer t_com;
    allreducePotential(prec, V_loc);
    t_com.stop();
    Printer::printTree(0, "Allreduce potential", this->getNNodes(NUMBER::Total), t_com.getWallTime());

    t_tot.stop();
    Printer::printFooter(0, t_tot, 2);
    Printer::setPrecision(oldprec);
}

void NuclearPotential::allreducePotential(double prec, QMFunction &V_loc) {
    QMFunction &V_tot = *this;

    // Add up local contributions into the grand master
    mpi::reduce_function(prec, V_loc, mpi::comm_orb);

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
}


/** @brief computes the interaction energy of the nuclear potential with a second set of nuclei.
 * 
 * @param[in] nucs the set of nuclei
 *
 * Note: this function is not suited to compute the nuclear self-energy
 *
 */
double NuclearOperator::trace(const Nuclei &nucs) {
    MSG_WARN("This routine has never been tested!");
    int nNucs = nucs.size();
    double E_nuc = 0.0;
    for (int i = 0; i < nNucs; i++) {
        const Nucleus &nuc_i = nucs[i];
        double Z_i = nuc_i.getCharge();
        const mrcpp::Coord<3> &R_i = nuc_i.getCoord();
        E_nuc += Z_i * this->r_m1.evalf(R_i);
    }
    return E_nuc;
}

} //namespace mrchem
