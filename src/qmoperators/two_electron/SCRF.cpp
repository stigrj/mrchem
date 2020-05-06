#include "SCRF.h"
#include "MRCPP/MWOperators"
#include "ReactionPotential.h"
#include "chemistry/Nucleus.h"
#include "chemistry/Permittivity.h"
#include "chemistry/chemistry_utils.h"
#include "qmfunctions/Density.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/QMFunction.h"
#include "qmfunctions/density_utils.h"
#include "qmfunctions/qmfunction_fwd.h"
#include "qmfunctions/qmfunction_utils.h"

using PoissonOperator_p = std::shared_ptr<mrcpp::PoissonOperator>;
using DerivativeOperator_p = std::shared_ptr<mrcpp::DerivativeOperator<3>>;
using OrbitalVector_p = std::shared_ptr<mrchem::OrbitalVector>;

namespace mrchem {
  SCRF::SCRF(Nuclei N, Permittivity e, OrbitalVector_p phi, PoissonOperator_p P, DerivativeOperator_p D, double orb_prec)
        : apply_prec(orb_prec)
        , epsilon(e)
        , poisson(P)
        , derivative(D)
        , rho_nuc(false)
        , rho_tot(false)
        , difference_potential(false)
        , potential(false) {
        rho_nuc = chemistry::compute_nuclear_density(this->apply_prec, N, 1000);
        updateTotalDensity(*phi, this->apply_prec);

        mrcpp::FunctionTree<3> *dx_cavity = new mrcpp::FunctionTree<3>(*MRA);
        mrcpp::FunctionTree<3> *dy_cavity = new mrcpp::FunctionTree<3>(*MRA);
        mrcpp::FunctionTree<3> *dz_cavity = new mrcpp::FunctionTree<3>(*MRA);
        d_cavity.push_back(std::make_tuple(1.0, dx_cavity));
        d_cavity.push_back(std::make_tuple(1.0, dy_cavity));
        d_cavity.push_back(std::make_tuple(1.0, dz_cavity));
        mrcpp::project<3>(this->apply_prec / 100, this->d_cavity, this->epsilon.getGradVector());
  }

  void SCRF::updateTotalDensity(OrbitalVector Phi,
                                double prec) { // pass the electron orbitals and computes the total density
    resetQMFunction(this->rho_tot);
    Density rho_el(false);
    density::compute(prec, rho_el, Phi, DensityType::Total);
    rho_el.rescale(-1.0);
    qmfunction::add(this->rho_tot, 1.0, rho_el, 1.0, this->rho_nuc, -1.0); // probably change this into a vector
}

QMFunction SCRF::updateGamma(QMFunction potential_nm1, double prec) {
    QMFunction gamma;
    gamma.alloc(NUMBER::Real);
    auto d_V = mrcpp::gradient(*derivative, potential_nm1.real());
    mrcpp::dot(prec, gamma.real(), d_V, d_cavity);
    gamma.rescale(std::log((epsilon.eps_in / epsilon.eps_out)) * (1.0 / (4.0 * MATHCONST::pi)));
    return gamma;
}

QMFunctionVector SCRF::makeTerms(double prec) {
    std::cout << "this potential:\t" << this->potential.integrate() << "\n";
    std::cout << "this difference_potential:\t" << this->difference_potential.integrate() << "\n";
    QMFunction vacuum_potential;
    QMFunction rho_eff;
    QMFunction eps;
    QMFunction eps_inv;
    QMFunction first_term;
    QMFunction gamma;
    QMFunction total_potential;

    vacuum_potential.alloc(NUMBER::Real);
    rho_eff.alloc(NUMBER::Real);
    eps.alloc(NUMBER::Real);
    eps_inv.alloc(NUMBER::Real);
    first_term.alloc(NUMBER::Real);
    total_potential.alloc(NUMBER::Real);
    mrcpp::apply(prec, vacuum_potential.real(), *poisson, this->rho_tot.real());

    epsilon.flipFunction(false);
    qmfunction::project(eps, epsilon, NUMBER::Real, prec / 100);
    epsilon.flipFunction(true);
    qmfunction::project(eps_inv, epsilon, NUMBER::Real, prec / 100);
    qmfunction::multiply(first_term, rho_tot, eps_inv, prec);
    qmfunction::add(rho_eff, 1.0, first_term, -1.0, rho_tot, -1.0);

    if (not this->potential.hasReal()) {
        QMFunction poisson_func;
        QMFunction V_n;
        poisson_func.alloc(NUMBER::Real);
        V_n.alloc(NUMBER::Real);

        gamma = updateGamma(vacuum_potential, prec);
        qmfunction::add(poisson_func, 1.0, rho_eff, 1.0, gamma, -1.0);
        mrcpp::apply(prec, V_n.real(), *poisson, poisson_func.real());
        this->difference_potential = V_n;
    }
    QMFunction Vr_np1;
    Vr_np1.alloc(NUMBER::Real);
    qmfunction::add(Vr_np1, 1.0, potential, 1.0, difference_potential, -1.0);

    qmfunction::add(total_potential, 1.0, Vr_np1, 1.0, vacuum_potential, -1.0);
    resetQMFunction(gamma);
    gamma = updateGamma(total_potential, prec);

    QMFunctionVector terms_vector;
    terms_vector.push_back(Vr_np1);
    terms_vector.push_back(gamma);
    terms_vector.push_back(rho_eff);
    terms_vector.push_back(vacuum_potential);

    return terms_vector;
}

void SCRF::updatePotential(QMFunction new_potential) {
    resetQMFunction(potential);
    qmfunction::deep_copy(this->potential, new_potential);
}

void SCRF::updateDifferencePotential(QMFunction diff_potential) {
    resetQMFunction(difference_potential);
    qmfunction::deep_copy(this->difference_potential, diff_potential);
    QMFunction new_potential;
    new_potential.alloc(NUMBER::Real);
    qmfunction::add(new_potential, 1.0, this->potential, 1.0, this->difference_potential, -1.0);
    this->Rp->updatePotential(new_potential);
}

double SCRF::getNuclearEnergy() {
    QMFunction V_n;
    QMFunction integral_product;
    V_n.alloc(NUMBER::Real);
    integral_product.alloc(NUMBER::Real);

    qmfunction::add(V_n, 1.0, this->potential, 1.0, this->difference_potential, -1.0);
    qmfunction::multiply(integral_product, this->rho_nuc, V_n, this->apply_prec);
    return integral_product.integrate().real();
}

double SCRF::getElectronicEnergy() {
    QMFunction V_n;
    QMFunction integral_product;
    QMFunction rho_el;
    V_n.alloc(NUMBER::Real);
    integral_product.alloc(NUMBER::Real);
    rho_el.alloc(NUMBER::Real);
    qmfunction::add(rho_el, 1.0, rho_tot, -1.0, rho_nuc, -1.0);
    qmfunction::add(V_n, 1.0, this->potential, 1.0, this->difference_potential, -1.0);
    qmfunction::multiply(integral_product, rho_el, V_n, this->apply_prec);
    return integral_product.integrate().real();
}

double SCRF::getTotalEnergy() {
    QMFunction V_n;
    QMFunction integral_product;
    V_n.alloc(NUMBER::Real);
    integral_product.alloc(NUMBER::Real);

    qmfunction::add(V_n, 1.0, this->potential, 1.0, this->difference_potential, -1.0);
    qmfunction::multiply(integral_product, this->rho_tot, V_n, this->apply_prec);
    return integral_product.integrate().real();
}

void SCRF::resetQMFunction(QMFunction &function) {
    if (function.hasReal()) function.free(NUMBER::Real);
    if (function.hasImag()) function.free(NUMBER::Imag);
    function.alloc(NUMBER::Real);
}
void SCRF::clear() {
  this->rho_tot.free(NUMBER::Real);
}
} // namespace mrchem