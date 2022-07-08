/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
 *
 * This file is part of MRChem.
 *
 * MRChem is free software: you can redistribute it and/or modify
 * it under the terms of the GNU Lesser General Public License as published by
 * the Free Software Foundation, either version 3 of the License, or
 * (at your option) any later version.
 *
 * MRChem is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public License
 * along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
 *
 * For information on the complete list of contributors to MRChem, see:
 * <https://mrchem.readthedocs.io/>
 */

#include "FockBuilder.h"

#include <MRCPP/Printer>
#include <MRCPP/Timer>

#include "CoulombOperator.h"
#include "ExchangeOperator.h"
#include "ReactionOperator.h"
#include "XCOperator.h"
#include "chemistry/chemistry_utils.h"
#include "properties/SCFEnergy.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/one_electron/ElectricFieldOperator.h"
#include "qmoperators/one_electron/IdentityOperator.h"
#include "qmoperators/one_electron/KineticOperator.h"
#include "qmoperators/one_electron/NablaOperator.h"
#include "qmoperators/one_electron/NuclearOperator.h"
#include "qmoperators/one_electron/ZoraOperator.h"
#include "qmoperators/qmoperator_utils.h"
#include "utils/math_utils.h"

using mrcpp::Printer;
using mrcpp::Timer;

namespace mrchem {

/** @brief build the Fock operator once all contributions are in place
 *
 */
void FockBuilder::build(double exx) {
    this->exact_exchange = exx;

    this->V = RankZeroOperator();
    if (this->nuc != nullptr) this->V += (*this->nuc);
    if (this->coul != nullptr) this->V += (*this->coul);
    if (this->ex != nullptr) this->V -= this->exact_exchange * (*this->ex);
    if (this->xc != nullptr) this->V += (*this->xc);
    if (this->ext != nullptr) this->V += (*this->ext);
    if (this->Ro != nullptr) this->V -= (*this->Ro);
}

/** @brief prepare operator for application
 *
 * @param prec: apply precision
 *
 * This will call the setup function of all underlying operators, and in particular
 * it will compute the internal exchange if there is an ExchangeOperator.
 */
void FockBuilder::setup(double prec) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    if (plevel == 2) {
        mrcpp::print::header(2, "Building Fock operator");
        mrcpp::print::value(2, "Precision", prec, "(rel)", 5);
        mrcpp::print::separator(2, '-');
    }
    if (this->mom != nullptr) this->momentum().setup(prec);
    this->potential().setup(prec);
    this->perturbation().setup(prec);

    if (isZora()) {
        Timer t_zora;
        auto c = getLightSpeed();
        mrcpp::print::header(3, "Building ZORA operators");
        mrcpp::print::value(3, "Precision", prec, "(rel)", 5);
        mrcpp::print::value(3, "Light speed", c, "(au)", 5);
        mrcpp::print::separator(3, '-');
        auto vz = collectZoraBasePotential();
        this->kappa = std::make_shared<ZoraOperator>(*vz, c, prec, false);
        this->kappa_inv = std::make_shared<ZoraOperator>(*vz, c, prec, true);
        this->zora_base = RankZeroOperator(vz);
        this->kappa->setup(prec);
        this->kappa_inv->setup(prec);
        this->zora_base.setup(prec);
        mrcpp::print::footer(3, t_zora, 2);
    };

    t_tot.stop();
    if (plevel == 2) mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Building Fock operator", t_tot);
}

/** @brief clear operator after application
 *
 * This will call the clear function of all underlying operators, and bring them back
 * to the state after construction. The operator can now be reused after another setup.
 */
void FockBuilder::clear() {
    if (this->mom != nullptr) this->momentum().clear();
    this->potential().clear();
    this->perturbation().clear();
    if (isZora()) {
        this->kappa->clear();
        this->kappa_inv->clear();
        this->zora_base.clear();
    }
}

/** @brief rotate orbitals of two-electron operators
 *
 * @param U: unitary transformation matrix
 *
 * This function should be used in case the orbitals are rotated *after* the FockBuilder
 * has been setup. In particular the ExchangeOperator needs to rotate the precomputed
 * internal exchange potentials.
 */
void FockBuilder::rotate(const ComplexMatrix &U) {
    if (this->ex != nullptr) this->ex->rotate(U);
}

/** @brief compute the SCF energy
 *
 * @param Phi: orbitals
 * @param F: Fock matrix
 *
 * This function will compute the total energy for a given OrbitalVector and
 * the corresponding Fock matrix. Tracing the kinetic energy operator is avoided
 * by tracing the Fock matrix and subtracting all other contributions.
 */
SCFEnergy FockBuilder::trace(OrbitalVector &Phi, const Nuclei &nucs) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Computing molecular energy");

    double E_kin = 0.0;  // Kinetic energy
    double E_nn = 0.0;   // Nuclear repulsion
    double E_en = 0.0;   // Nuclear-electronic interaction
    double E_ee = 0.0;   // Electronic repulsion
    double E_x = 0.0;    // Exact Exchange
    double E_xc = 0.0;   // Exchange and Correlation
    double E_eext = 0.0; // External field contribution to the electronic energy
    double E_next = 0.0; // External field contribution to the nuclear energy
    double Er_nuc = 0.0; // Nuclear reaction energy
    double Er_el = 0.0;  // Electronic reaction energy
    double Er_tot = 0.0; // Total reaction energy

    // Nuclear part
    if (this->nuc != nullptr) E_nn = chemistry::compute_nuclear_repulsion(nucs);
    if (this->ext != nullptr) E_next = -this->ext->trace(nucs).real();

    // Reaction potential part
    if (this->Ro != nullptr) {
        Er_nuc = 0.5 * this->Ro->getNuclearEnergy();
        Er_tot = 0.5 * this->Ro->getTotalEnergy();
        Er_el = 0.5 * this->Ro->getElectronicEnergy();
    }

    // Kinetic part
    if (isZora()) {
        E_kin = qmoperator::calc_kinetic_trace(momentum(), *this->kappa, Phi).real();
    } else {
        E_kin = qmoperator::calc_kinetic_trace(momentum(), Phi);
    }

    // Electronic part
    if (this->nuc != nullptr) E_en = this->nuc->trace(Phi).real();
    if (this->coul != nullptr) E_ee = 0.5 * this->coul->trace(Phi).real();
    if (this->ex != nullptr) E_x = -this->exact_exchange * this->ex->trace(Phi).real();
    if (this->xc != nullptr) E_xc = this->xc->getEnergy();
    if (this->ext != nullptr) E_eext = this->ext->trace(Phi).real();
    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Computing molecular energy", t_tot);

    return SCFEnergy{E_kin, E_nn, E_en, E_ee, E_x, E_xc, E_next, E_eext, Er_tot, Er_nuc, Er_el};
}

ComplexMatrix FockBuilder::operator()(OrbitalVector &bra, OrbitalVector &ket) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Computing Fock matrix");

    ComplexMatrix T_mat = ComplexMatrix::Zero(bra.size(), ket.size());
    if (isZora()) {
        T_mat = qmoperator::calc_kinetic_matrix(momentum(), *this->kappa, bra, ket);
    } else {
        T_mat = qmoperator::calc_kinetic_matrix(momentum(), bra, ket);
    }

    ComplexMatrix V_mat = ComplexMatrix::Zero(bra.size(), ket.size());
    V_mat += potential()(bra, ket);

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Computing Fock matrix", t_tot);
    return T_mat + V_mat;
}

OrbitalVector FockBuilder::buildHelmholtzArgument(double prec, OrbitalVector Phi, ComplexMatrix F_mat, ComplexMatrix L_mat) {
    Timer t_tot;
    auto plevel = Printer::getPrintLevel();
    mrcpp::print::header(2, "Computing Helmholtz argument");

    Timer t_rot;
    OrbitalVector Psi = orbital::rotate(Phi, L_mat - F_mat, prec);
    mrcpp::print::time(2, "Rotating orbitals", t_rot);

    OrbitalVector out;
    if (isZora()) {
        out = buildHelmholtzArgumentZORA(Phi, Psi, F_mat.real().diagonal(), prec);
    } else {
        out = buildHelmholtzArgumentNREL(Phi, Psi);
    }
    Psi.clear();

    mrcpp::print::footer(2, t_tot, 2);
    if (plevel == 1) mrcpp::print::time(1, "Computing Helmholtz argument", t_tot);
    return out;
}

// Take 1 in notes on Overleaf
OrbitalVector FockBuilder::buildHelmholtzArgumentZORA(OrbitalVector &Phi, OrbitalVector &Psi, DoubleVector eps, double prec) {
    // Get necessary operators
    double c = getLightSpeed();
    double two_cc = 2.0 * c * c;
    MomentumOperator &p = momentum();
    RankZeroOperator &V = potential();
    RankZeroOperator &kappa = *this->kappa;
    RankZeroOperator &kappa_m1 = *this->kappa_inv;
    RankZeroOperator &V_zora = this->zora_base;

    RankZeroOperator operOne = 0.5 * tensor::dot(p(kappa), p);
    RankZeroOperator operThree = kappa * V_zora;
    operOne.setup(prec);
    operThree.setup(prec);

    // Compute OrbitalVectors
    Timer t_1;
    OrbitalVector termOne = operOne(Phi);
    mrcpp::print::time(2, "Computing gradient term", t_1);

    Timer t_2;
    OrbitalVector termTwo = V(Phi);
    mrcpp::print::time(2, "Computing potential term", t_2);

    // Compute transformed orbitals scaled by diagonal Fock elements
    Timer t_3;
    OrbitalVector epsPhi = orbital::deep_copy(Phi);
    for (int i = 0; i < epsPhi.size(); i++) {
        if (not mpi::my_orb(epsPhi[i])) continue;
        epsPhi[i].rescale(eps[i] / two_cc);
    }
    OrbitalVector termThree = operThree(epsPhi);
    mrcpp::print::time(2, "Computing rescaled potential term", t_3);

    auto normsOne = orbital::get_norms(termOne);
    auto normsTwo = orbital::get_norms(termTwo);
    auto normsThree = orbital::get_norms(termThree);
    auto normsPsi = orbital::get_norms(Psi);

    // Add up all the terms
    Timer t_add;
    OrbitalVector arg = orbital::deep_copy(termOne);
    for (int i = 0; i < arg.size(); i++) {
        if (not mpi::my_orb(arg[i])) continue;
        arg[i].add(1.0, termTwo[i]);
        arg[i].add(1.0, termThree[i]);
        arg[i].add(1.0, Psi[i]);
    }
    mrcpp::print::time(2, "Adding contributions", t_add);

    operThree.clear();
    operOne.clear();

    Timer t_kappa;
    auto out = kappa_m1(arg);
    mrcpp::print::time(2, "Applying kappa inverse", t_kappa);
    return out;
}

// Non-relativistic Helmholtz argument
OrbitalVector FockBuilder::buildHelmholtzArgumentNREL(OrbitalVector &Phi, OrbitalVector &Psi) {
    // Get necessary operators
    RankZeroOperator &V = this->potential();

    // Compute OrbitalVectors
    Timer t_pot;
    OrbitalVector termOne = V(Phi);
    mrcpp::print::time(2, "Computing potential term", t_pot);

    // Add up all the terms
    Timer t_add;
    OrbitalVector out = orbital::deep_copy(termOne);
    for (int i = 0; i < out.size(); i++) {
        if (not mpi::my_orb(out[i])) continue;
        out[i].add(1.0, Psi[i]);
    };
    mrcpp::print::time(2, "Adding contributions", t_add);
    return out;
}

void FockBuilder::setZoraType(bool has_nuc, bool has_coul, bool has_xc) {
    this->zora_has_nuc = has_nuc;
    this->zora_has_coul = has_coul;
    this->zora_has_xc = has_xc;
}

std::shared_ptr<QMPotential> FockBuilder::collectZoraBasePotential() {
    Timer timer;
    auto vz = std::make_shared<QMPotential>(1, false);
    if (zora_has_nuc) {
        if (getNuclearOperator() != nullptr) {
            auto &vnuc = static_cast<QMPotential &>(getNuclearOperator()->getRaw(0, 0));
            if (not vnuc.hasReal()) MSG_ERROR("ZORA: Adding empty nuclear potential");
            vz->add(1.0, vnuc);
        } else {
            MSG_ERROR("ZORA: Nuclear requested but not available");
        }
    }
    if (zora_has_coul) {
        if (getCoulombOperator() != nullptr) {
            auto &coul = static_cast<QMPotential &>(getCoulombOperator()->getRaw(0, 0));
            if (not coul.hasReal()) MSG_INFO("ZORA: Adding empty Coulomb potential");
            vz->add(1.0, coul);
        } else {
            MSG_ERROR("ZORA: Coulomb requested but not available");
        }
    }
    if (zora_has_xc) {
        if (getXCOperator() != nullptr) {
            getXCOperator()->setSpin(SPIN::Alpha);
            auto &xc = static_cast<QMPotential &>(getXCOperator()->getRaw(0, 0));
            if (not xc.hasReal()) MSG_ERROR("ZORA: Adding empty XC potential");
            vz->add(1.0, xc);
            getXCOperator()->clearSpin();
        } else {
            MSG_ERROR("ZORA: XC requested but not available");
        }
    }
    print_utils::qmfunction(2, "ZORA operator (base)", *vz, timer);
    return vz;
}

} // namespace mrchem
