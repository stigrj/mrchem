#include "MRCPP/Printer"
#include "MRCPP/Timer"

#include "parallel.h"
#include "utils/math_utils.h"
#include "utils/RRMaximizer.h"

#include "qmfunctions.h"
#include "Orbital.h"
#include "NuclearCorrelationOperator.h"

using mrcpp::Timer;
using mrcpp::Printer;
using mrcpp::FunctionTree;
using mrcpp::FunctionTreeVector;

namespace mrchem {
extern mrcpp::MultiResolutionAnalysis<3> *MRA; // Global MRA

namespace orbital {

/****************************************
 * Orbital related standalone functions *
 ****************************************/

/** @brief Compute <bra|ket> = int bra^\dag(r) * ket(r) dr.
 *
 *  Complicated by the fact that both bra and ket can be interpreted
 *  as complex conjugate versions of themselves. Notice that the <bra|
 *  position is already complex conjugated.
 *
 */
ComplexDouble dot(Orbital bra, Orbital ket, NuclearCorrelationOperator *R) {
    if ((bra.spin() == SPIN::Alpha) and (ket.spin() == SPIN::Beta)) return 0.0;
    if ((bra.spin() == SPIN::Beta) and (ket.spin() == SPIN::Alpha)) return 0.0;


    Orbital RRket = ket;
    if (R != nullptr) RRket = R->getR2()(ket);

    double bra_conj(1.0), ket_conj(1.0);
    if (bra.conjugate()) bra_conj = -1.0;
    if (ket.conjugate()) ket_conj = -1.0;

    double rr(0.0), ri(0.0), ir(0.0), ii(0.0);
    if (bra.hasReal() and ket.hasReal()) rr = mrcpp::dot(bra.real(), RRket.real());
    if (bra.hasReal() and ket.hasImag()) ri = mrcpp::dot(bra.real(), RRket.imag());
    if (bra.hasImag() and ket.hasReal()) ir = mrcpp::dot(bra.imag(), RRket.real());
    if (bra.hasImag() and ket.hasImag()) ii = mrcpp::dot(bra.imag(), RRket.imag());
    if (R != nullptr) RRket.free();

    double real_part = rr + bra_conj*ket_conj*ii;
    double imag_part = ket_conj*ri - bra_conj*ir;
    return ComplexDouble(real_part, imag_part);
}

/** @brief Compare spin and occupancy of two orbitals
 *
 *  Returns true if orbital parameters are the same.
 *
 */
bool compare(const Orbital &orb_a, const Orbital &orb_b) {
    bool comp = true;
    if (compare_occ(orb_a, orb_b) < 0) {
        MSG_WARN("Different occupancy");
        comp = false;
    }
    if (compare_spin(orb_a, orb_b) < 0) {
        MSG_WARN("Different spin");
        comp = false;
    }
    return comp;
}

/** @brief Compare occupancy of two orbitals
 *
 *  Returns the common occupancy if they match, -1 if they differ.
 *
 */
int compare_occ(const Orbital &orb_a, const Orbital &orb_b) {
    int comp = -1;
    if (orb_a.occ() == orb_b.occ()) comp = orb_a.occ();
    return comp;
}

/** @brief Compare spin of two orbitals
 *
 *  Returns the common spin if they match, -1 if they differ.
 *
 */
int compare_spin(const Orbital &orb_a, const Orbital &orb_b) {
    int comp = -1;
    if (orb_a.spin() == orb_b.spin()) comp = orb_a.spin();
    return comp;
}

/** @brief out = a*inp_a + b*inp_b
  * Complicated by the fact that both inputs can be interpreted as complex
  * conjugate versions of themselves. */
Orbital add(ComplexDouble a, Orbital inp_a, ComplexDouble b, Orbital inp_b, double prec) {
    ComplexVector coefs(2);
    coefs(0) = a;
    coefs(1) = b;

    OrbitalVector orbs;
    orbs.push_back(inp_a);
    orbs.push_back(inp_b);

    return multiply(coefs, orbs, prec);
}

/** @brief out_i = a*(inp_a)_i + b*(inp_b)_i
 *
 *  Component-wise addition of orbitals.
 *
 */
OrbitalVector add(ComplexDouble a, OrbitalVector &inp_a,
                  ComplexDouble b, OrbitalVector &inp_b,
                  double prec) {
    if (inp_a.size() != inp_b.size()) MSG_ERROR("Size mismatch");

    OrbitalVector out;
    for (int i = 0; i < inp_a.size(); i++) {
        Orbital out_i = add(a, inp_a[i], b, inp_b[i]);
        out.push_back(out_i);
    }
    return out;
}

/** @brief out = inp_a * inp_b
  *
  * Complicated by the fact that both inputs can be interpreted
  * as complex conjugate versions of themselves.
  *
  */
Orbital multiply(Orbital inp_a, Orbital inp_b, double prec) {
    int occ = compare_occ(inp_a, inp_b);
    int spin = compare_spin(inp_a, inp_b);
    Orbital out(spin, occ);

    double a_conj(1.0), b_conj(1.0);
    if (inp_a.conjugate()) a_conj = -1.0;
    if (inp_b.conjugate()) b_conj = -1.0;

    { // Real part
        FunctionTreeVector<3> vec;
        if (inp_a.hasReal() and inp_b.hasReal()) {
            FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
            double coef = 1.0;
            if (prec < 0.0) {
                // Union grid
                mrcpp::build_grid(*tree, inp_a.real());
                mrcpp::build_grid(*tree, inp_b.real());
                mrcpp::multiply(prec, *tree, coef, inp_a.real(), inp_b.real(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, inp_a.real(), inp_b.real());
            }
            vec.push_back(std::make_tuple(1.0, tree));
        }
        if (inp_a.hasImag() and inp_b.hasImag()) {
            FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
            double coef = -1.0*a_conj*b_conj;
            if (prec < 0.0) {
                mrcpp::build_grid(*tree, inp_a.imag());
                mrcpp::build_grid(*tree, inp_b.imag());
                mrcpp::multiply(prec, *tree, coef, inp_a.imag(), inp_b.imag(), 0);
            } else {
                mrcpp::multiply(prec, *tree, coef, inp_a.imag(), inp_b.imag());
            }
            vec.push_back(std::make_tuple(1.0, tree));
        }
        if (vec.size() == 1) {
            out.setReal(&mrcpp::get_func(vec, 0));
            mrcpp::clear(vec, false);
        }
        if (vec.size() == 2) {
            out.alloc(NUMBER::Real);
            mrcpp::build_grid(out.real(), vec);
            mrcpp::add(prec, out.real(), vec, 0);
            mrcpp::clear(vec, true);
        }
    }

    { // Imaginary part
        FunctionTreeVector<3> vec;
        if (inp_a.hasReal() and inp_b.hasImag()) {
            FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
            double coef = b_conj;
            if (prec < 0.0) {
                // Union grid
                mrcpp::build_grid(*tree, inp_a.real());
                mrcpp::build_grid(*tree, inp_b.imag());
                mrcpp::multiply(prec, *tree, coef, inp_a.real(), inp_b.imag(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, inp_a.real(), inp_b.imag());
            }
            vec.push_back(std::make_tuple(1.0, tree));
        }
        if (inp_a.hasImag() and inp_b.hasReal()) {
            FunctionTree<3> *tree = new FunctionTree<3>(*MRA);
            double coef = a_conj;
            if (prec < 0.0) {
                // Union grid
                mrcpp::build_grid(*tree, inp_a.imag());
                mrcpp::build_grid(*tree, inp_b.real());
                mrcpp::multiply(prec, *tree, coef, inp_a.imag(), inp_b.real(), 0);
            } else {
                // Adaptive grid
                mrcpp::multiply(prec, *tree, coef, inp_a.imag(), inp_b.real());
            }
            vec.push_back(std::make_tuple(1.0, tree));
        }
        if (vec.size() == 1) {
            out.setImag(&mrcpp::get_func(vec, 0));
            mrcpp::clear(vec, false);
        }
        if (vec.size() == 2) {
            out.alloc(NUMBER::Imag);
            mrcpp::build_grid(out.imag(), vec);
            mrcpp::add(prec, out.imag(), vec, 0);
            mrcpp::clear(vec, true);
        }
    }

    return out;
}

/** @brief out = c_0*inp_0 + c_1*inp_1 + ...
  *
  * Complicated by the fact that both inputs can be interpreted as complex
  * conjugate versions of themselves.
  *
  */
Orbital multiply(const ComplexVector &c, OrbitalVector &inp, double prec) {
    if (c.size() != inp.size()) MSG_ERROR("Size mismatch");
    double thrs = mrcpp::MachineZero;

    Orbital out;
    // set output spin from first contributing input
    for (int i = 0; i < inp.size(); i++) {
        if (std::abs(c[i]) < thrs) continue;
        out = inp[i].paramCopy();
        break;
    }
    // all contributing input spins must be equal
    for (int i = 0; i < inp.size(); i++) {
        if (std::abs(c[i]) < thrs) continue;
        if (out.spin() != inp[i].spin()) {
            // empty orbitals with wrong spin can occur
            if (inp[i].hasReal()) MSG_FATAL("Mixing spins");
            if (inp[i].hasImag()) MSG_FATAL("Mixing spins");
        }
    }

    FunctionTreeVector<3> rvec;
    FunctionTreeVector<3> ivec;

    for (int i = 0; i < inp.size(); i++) {
        bool cHasReal = (std::abs(c[i].real()) > thrs);
        bool cHasImag = (std::abs(c[i].imag()) > thrs);

        double conj(1.0);
        if (inp[i].conjugate()) conj = -1.0;

        if (cHasReal and inp[i].hasReal()) rvec.push_back(std::make_tuple(      c[i].real(), &inp[i].real()));
        if (cHasImag and inp[i].hasImag()) rvec.push_back(std::make_tuple(-conj*c[i].imag(), &inp[i].imag()));

        if (cHasImag and inp[i].hasReal()) ivec.push_back(std::make_tuple(      c[i].imag(), &inp[i].real()));
        if (cHasReal and inp[i].hasImag()) ivec.push_back(std::make_tuple( conj*c[i].real(), &inp[i].imag()));
    }

    if (rvec.size() > 0) {
        out.alloc(NUMBER::Real);
        if (prec < 0.0) {
            mrcpp::build_grid(out.real(), rvec);
            mrcpp::add(prec, out.real(), rvec, 0);
        } else {
            mrcpp::add(prec, out.real(), rvec);
        }
    }
    if (ivec.size() > 0) {
        out.alloc(NUMBER::Imag);
        if (prec < 0.0) {
            mrcpp::build_grid(out.imag(), ivec);
            mrcpp::add(prec, out.imag(), ivec, 0);
        } else {
            mrcpp::add(prec, out.imag(), ivec);
        }
    }
    return out;
}

/** @brief Orbital transformation out_vec = U*inp_vec
 *
 * The transformation matrix is not necessarily square.
 *
 */
OrbitalVector multiply(const ComplexMatrix &U, OrbitalVector &inp, double prec) {
    if (mpi::orb_size > 1) NOT_IMPLEMENTED_ABORT;

    OrbitalVector out;
    for (int i = 0; i < U.rows(); i++) {
        const ComplexVector &c = U.row(i);
        Orbital out_i = multiply(c, inp, prec);
        out.push_back(out_i);
    }
    return out;
}

/** @brief Deep copy
 *
 * New orbitals are constructed as deep copies of the input set.
 *
 */
OrbitalVector deep_copy(OrbitalVector &inp) {
    OrbitalVector out;
    for (int i = 0; i < inp.size(); i++) {
        Orbital out_i = inp[i].deepCopy();
        out.push_back(out_i);
    }
    return out;
}

/** @brief Parameter copy
 *
 * New orbitals are constructed as parameter copies of the input set.
 *
 */
OrbitalVector param_copy(const OrbitalVector &inp) {
    OrbitalVector out;
    for (int i = 0; i < inp.size(); i++) {
        Orbital out_i = inp[i].paramCopy();
        out.push_back(out_i);
    }
    return out;
}

/** @brief Adjoin two vectors
 *
 * The orbitals of the input vector are appended to
 * (*this) vector, the ownership is transferred. Leaves
 * the input vector empty.
 *
 */
OrbitalVector adjoin(OrbitalVector &inp_a, OrbitalVector &inp_b) {
    OrbitalVector out;
    for (int i = 0; i < inp_a.size(); i++) out.push_back(inp_a[i]);
    for (int i = 0; i < inp_b.size(); i++) out.push_back(inp_b[i]);
    inp_a.clear();
    inp_b.clear();
    return out;
}

/** @brief Disjoin vector in two parts
 *
 * All orbitals of a particular spin is collected in a new vector
 * and returned. These orbitals are removed from (*this) vector,
 * and the ownership is transferred.
 *
 */
OrbitalVector disjoin(OrbitalVector &inp, int spin) {
    OrbitalVector out;
    OrbitalVector tmp;
    for (int i = 0; i < inp.size(); i++) {
        Orbital &orb_i = inp[i];
        if (orb_i.spin() == spin) {
            out.push_back(orb_i);
        } else {
            tmp.push_back(orb_i);
        }
    }
    inp.clear();
    inp = tmp;
    return out;
}

/** @brief Write orbitals to disk
 *
 * @param Phi: orbitals to save
 * @param file: file name prefix
 * @param n_orbs: number of orbitals to save
 *
 * The given file name (e.g. "phi") will be appended with orbital number ("phi_0").
 * Produces separate files for meta data ("phi_0.meta"), real ("phi_0_re.tree") and
 * imaginary ("phi_0_im.tree") parts. Negative n_orbs means that all orbitals in the
 * vector are saved.
 */
void save_orbitals(OrbitalVector &Phi, const std::string &file, int n_orbs) {
    if (n_orbs < 0) n_orbs = Phi.size();
    if (n_orbs > Phi.size()) MSG_ERROR("Index out of bounds");
    for (int i = 0; i < n_orbs; i++) {
        std::stringstream orbname;
        orbname << file << "_" << i;
        Phi[i].saveOrbital(orbname.str());
    }
}

/** @brief Read orbitals from disk
 *
 * @param file: file name prefix
 * @param n_orbs: number of orbitals to read
 *
 * The given file name (e.g. "phi") will be appended with orbital number ("phi_0").
 * Reads separate files for meta data ("phi_0.meta"), real ("phi_0_re.tree") and
 * imaginary ("phi_0_im.tree") parts. Negative n_orbs means that all orbitals matching
 * the prefix name will be read.
 */
OrbitalVector load_orbitals(const std::string &file, int n_orbs) {
    OrbitalVector Phi;
    for (int i = 0; true; i++) {
        if (n_orbs > 0 and i >= n_orbs) break;
        std::stringstream orbname;
        orbname << file << "_" << i;
        Orbital phi_i;
        phi_i.loadOrbital(orbname.str());
        if (phi_i.hasReal() or phi_i.hasImag()) {
            Phi.push_back(phi_i);
        } else {
            break;
        }
    }
    return Phi;
}

/** @brief Frees each orbital in the vector
 *
 * Leaves an empty vector. Orbitals are freed.
 *
 */
void free(OrbitalVector &vec) {
    for (int i = 0; i < vec.size(); i++) vec[i].free();
    vec.clear();
}

/** @brief Normalize all orbitals in the set */
void normalize(OrbitalVector &vec) {
    for (int i = 0; i < vec.size(); i++) {
        vec[i].normalize();
    }
}

/** @brief Gram-Schmidt orthogonalize orbitals within the set */
void orthogonalize(OrbitalVector &vec) {
    for (int i = 0; i < vec.size(); i++) {
        Orbital &orb_i = vec[i];
        for (int j = 0; j < i; j++) {
            Orbital &orb_j = vec[j];
            orb_i.orthogonalize(orb_j);
        }
    }
}

/** @brief Orthogonalize the out orbital against all orbitals in inp */
void orthogonalize(OrbitalVector &vec, OrbitalVector &inp) {
    for (int i = 0; i < vec.size(); i++) {
        vec[i].orthogonalize(inp);
    }
}

ComplexMatrix calc_overlap_matrix(OrbitalVector &braket, NuclearCorrelationOperator *R) {
    return orbital::calc_overlap_matrix(braket, braket, R);
}

ComplexMatrix calc_overlap_matrix(OrbitalVector &bra,
                                  OrbitalVector &ket,
                                  NuclearCorrelationOperator *R) {
    int Ni = bra.size();
    int Nj = ket.size();
    ComplexMatrix S(Ni, Nj);
    S.setZero();

    for (int i = 0; i < Ni; i++) {
        for (int j = 0; j < Nj; j++) {
            S(i,j) = orbital::dot(bra[i], ket[j], R);
        }
    }
    return S;
}

/** @brief Compute Löwdin orthonormalization matrix
 *
 * @param Phi: orbitals to orthonomalize
 *
 * Computes the inverse square root of the orbital overlap matrix S^(-1/2)
 */
ComplexMatrix calc_lowdin_matrix(OrbitalVector &Phi, NuclearCorrelationOperator *R) {
    Timer timer;
    printout(1, "Calculating Löwdin orthonormalization matrix      ");

    ComplexMatrix S_tilde = orbital::calc_overlap_matrix(Phi, R);
    ComplexMatrix S_m12 = math_utils::hermitian_matrix_pow(S_tilde, -1.0/2.0);

    timer.stop();
    println(1, timer.getWallTime());
    return S_m12;
}

/** @brief Minimize the spatial extension of orbitals, by orbital rotation
 *
 * @param Phi: orbitals to localize
 *
 * Minimizes \f$  \sum_{i=1,N}\langle i| {\bf R^2}  | i \rangle - \langle i| {\bf R}| i \rangle^2 \f$
 * which is equivalent to maximizing \f$  \sum_{i=1,N}\langle i| {\bf R}| i \rangle^2\f$
 *
 * The resulting transformation includes the orthonormalization of the orbitals.
 * Orbitals are rotated in place, and the transformation matrix is returned.
 */
ComplexMatrix localize(double prec,
                       OrbitalVector &Phi,
                       NuclearCorrelationOperator *R) {
    Printer::printHeader(0, "Localizing orbitals");
    Timer timer;

    ComplexMatrix U;
    int n_it = 0;
    if (Phi.size() > 1) {
        Timer rr_t;
        RRMaximizer rr(prec, Phi, R);
        n_it = rr.maximize();
        rr_t.stop();
        Printer::printDouble(0, "Computing Foster-Boys matrix", rr_t.getWallTime(), 5);

        if (n_it > 0) {
            println(0, " Converged after iteration   " << std::setw(30) << n_it);
            U = rr.getTotalU().transpose().cast<ComplexDouble>();
        } else {
            println(0, " Foster-Boys localization did not converge!");
        }
    } else {
        println(0, " Cannot localize less than two orbitals");
    }

    if (n_it <= 0) {
        Timer orth_t;
        U = orbital::calc_lowdin_matrix(Phi, R);
        orth_t.stop();
        Printer::printDouble(0, "Computing Lowdin matrix", orth_t.getWallTime(), 5);
    }

    Timer rot_t;
    OrbitalVector Psi = orbital::multiply(U, Phi, prec);
    orbital::free(Phi);
    Phi = Psi;
    rot_t.stop();
    Printer::printDouble(0, "Rotating orbitals", rot_t.getWallTime(), 5);

    timer.stop();
    Printer::printFooter(0, timer, 2);
    return U;
}

/** @brief Perform the orbital rotation that diagonalizes the Fock matrix
 *
 * @param Phi: orbitals to rotate
 * @param F: Fock matrix to diagonalize
 *
 * The resulting transformation includes the orthonormalization of the orbitals.
 * Orbitals are rotated in place and Fock matrix is diagonalized in place.
 * The transformation matrix is returned.
 */
ComplexMatrix diagonalize(double prec,
                          OrbitalVector &Phi,
                          ComplexMatrix &F,
                          NuclearCorrelationOperator *R) {
    Printer::printHeader(0, "Digonalizing Fock matrix");
    Timer timer;

    Timer orth_t;
    ComplexMatrix S_m12 = orbital::calc_lowdin_matrix(Phi, R);
    F = S_m12.transpose()*F*S_m12;
    orth_t.stop();
    Printer::printDouble(0, "Computing Lowdin matrix", orth_t.getWallTime(), 5);

    Timer diag_t;
    ComplexMatrix U = ComplexMatrix::Zero(F.rows(), F.cols());
    int np = orbital::size_paired(Phi);
    int na = orbital::size_alpha(Phi);
    int nb = orbital::size_beta(Phi);
    if (np > 0) math_utils::diagonalize_block(F, U, 0,       np);
    if (na > 0) math_utils::diagonalize_block(F, U, np,      na);
    if (nb > 0) math_utils::diagonalize_block(F, U, np + na, nb);
    U = U * S_m12;
    diag_t.stop();
    Printer::printDouble(0, "Diagonalizing matrix", diag_t.getWallTime(), 5);

    Timer rot_t;
    OrbitalVector Psi = orbital::multiply(U, Phi, prec);
    orbital::free(Phi);
    Phi = Psi;
    rot_t.stop();
    Printer::printDouble(0, "Rotating orbitals", rot_t.getWallTime(), 5);

    timer.stop();
    Printer::printFooter(0, timer, 2);
    return U;
}

/** @brief Perform the Löwdin orthonormalization
 *
 * @param Phi: orbitals to orthonormalize
 *
 * Orthonormalizes the orbitals by multiplication of the Löwdin matrix S^(-1/2).
 * Orbitals are rotated in place, and the transformation matrix is returned.
 */
ComplexMatrix orthonormalize(double prec,
                             OrbitalVector &Phi,
                             NuclearCorrelationOperator *R) {
    ComplexMatrix U = orbital::calc_lowdin_matrix(Phi, R);
    OrbitalVector Psi = orbital::multiply(U, Phi, prec);
    orbital::free(Phi);
    Phi = Psi;
    return U;
}

/** @brief Returns the number of occupied orbitals */
int size_occupied(const OrbitalVector &vec) {
    int nOcc = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].occ() > 0) nOcc++;
    }
    return nOcc;
}

/** @brief Returns the number of empty orbitals */
int size_empty(const OrbitalVector &vec) {
    int nEmpty = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].occ() == 0) nEmpty++;
    }
    return nEmpty;
}

/** @brief Returns the number of singly occupied orbitals */
int size_singly(const OrbitalVector &vec) {
    int nSingly = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].occ() == 1) nSingly++;
    }
    return nSingly;
}

/** @brief Returns the number of doubly occupied orbitals */
int size_doubly(const OrbitalVector &vec) {
    int nDoubly = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].occ() == 1) nDoubly++;
    }
    return nDoubly;
}

/** @brief Returns the number of paired orbitals */
int size_paired(const OrbitalVector &vec) {
    int nPaired = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].spin() == SPIN::Paired) nPaired++;
    }
    return nPaired;
}

/** @brief Returns the number of alpha orbitals */
int size_alpha(const OrbitalVector &vec) {
    int nAlpha = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].spin() == SPIN::Alpha) nAlpha++;
    }
    return nAlpha;
}

/** @brief Returns the number of beta orbitals */
int size_beta(const OrbitalVector &vec) {
    int nBeta = 0;
    for (int i = 0; i < vec.size(); i++) {
        if (vec[i].spin() == SPIN::Beta) nBeta++;
    }
    return nBeta;
}

/** @brief Returns the spin multiplicity of the vector */
int get_multiplicity(const OrbitalVector &vec) {
    int nAlpha = get_electron_number(vec, SPIN::Alpha);
    int nBeta = get_electron_number(vec, SPIN::Beta);
    int S = std::abs(nAlpha - nBeta);
    return S + 1;
}

/** @brief Returns the number of electrons with the given spin
 *
 * Paired spin (default input) returns the total number of electrons.
 *
 */
int get_electron_number(const OrbitalVector &vec, int spin) {
    int nElectrons = 0;
    for (int i = 0; i < vec.size(); i++) {
        const Orbital &orb = vec[i];
        if (spin == SPIN::Paired) {
            nElectrons += orb.occ();
        } else if (spin == SPIN::Alpha) {
            if (orb.spin() == SPIN::Paired or orb.spin() == SPIN::Alpha) {
                nElectrons += 1;
            }
        } else if (spin == SPIN::Beta) {
            if (orb.spin() == SPIN::Paired or orb.spin() == SPIN::Beta) {
                nElectrons += 1;
            }
        } else {
            MSG_ERROR("Invalid spin argument");
        }
    }
    return nElectrons;
}

/** @brief Returns a vector containing the orbital errors */
DoubleVector get_errors(const OrbitalVector &vec) {
    int nOrbs = vec.size();
    DoubleVector errors = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        errors(i) = vec[i].error();
    }
    return errors;
}

/** @brief Assign errors to each orbital.
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void set_errors(OrbitalVector &vec, const DoubleVector &errors) {
    if (vec.size() != errors.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < vec.size(); i++) {
        vec[i].setError(errors(i));
    }
}

/** @brief Returns a vector containing the orbital spins */
IntVector get_spins(const OrbitalVector &vec) {
    int nOrbs = vec.size();
    IntVector spins = IntVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        spins(i) = vec[i].spin();
    }
    return spins;
}

/** @brief Assigns spin to each orbital
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void set_spins(OrbitalVector &vec, const IntVector &spins) {
    if (vec.size() != spins.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < vec.size(); i++) {
        vec[i].setSpin(spins(i));
    }
}

/** @brief Returns a vector containing the orbital occupancies */
IntVector get_occupancies(const OrbitalVector &vec) {
    int nOrbs = vec.size();
    IntVector occ = IntVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        occ(i) = vec[i].occ();
    }
    return occ;
}

/** @brief Assigns spin to each orbital
 *
 * Length of input vector must match the number of orbitals in the set.
 *
 */
void set_occupancies(OrbitalVector &vec, const IntVector &occ) {
    if (vec.size() != occ.size()) MSG_ERROR("Size mismatch");
    for (int i = 0; i < vec.size(); i++) {
        vec[i].setOcc(occ(i));
    }
}

/** @brief Returns a vector containing the orbital square norms */
DoubleVector get_squared_norms(const OrbitalVector &vec) {
    int nOrbs = vec.size();
    DoubleVector norms = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        if (mpi::my_orb(vec[i])) norms(i) = vec[i].squaredNorm();
    }
    return norms;
}

/** @brief Returns a vector containing the orbital norms */
DoubleVector get_norms(const OrbitalVector &vec) {
    int nOrbs = vec.size();
    DoubleVector norms = DoubleVector::Zero(nOrbs);
    for (int i = 0; i < nOrbs; i++) {
        if (mpi::my_orb(vec[i])) norms(i) = vec[i].norm();
    }
    return norms;
}

void print(const OrbitalVector &vec) {
    Printer::setScientific();
    printout(0, "============================================================\n");
    printout(0, " OrbitalVector:");
    printout(0, std::setw(4) << vec.size()          << " orbitals  ");
    printout(0, std::setw(4) << size_occupied(vec)  << " occupied  ");
    printout(0, std::setw(4) << get_electron_number(vec) << " electrons\n");
    printout(0, "------------------------------------------------------------\n");
    printout(0, "   n   Rank        Norm          Nodes  Spin Occ     Error  \n");
    printout(0, "------------------------------------------------------------\n");
    for (int i = 0; i < vec.size(); i++) {
        println(0, std::setw(4) << i << vec[i]);
    }
    printout(0, "============================================================\n\n\n");
}

} //namespace orbital


/****************************************
 * Density related standalone functions *
 ****************************************/

void density::compute(double prec,
                      Density &rho,
                      Orbital phi,
                      int spin,
                      NuclearCorrelationOperator *R) {
    double occ_a(0.0), occ_b(0.0), occ_p(0.0);
    if (phi.spin() == SPIN::Alpha)  occ_a = (double) phi.occ();
    if (phi.spin() == SPIN::Beta)   occ_b = (double) phi.occ();
    if (phi.spin() == SPIN::Paired) occ_p = (double) phi.occ();

    double occ(0.0);
    if (spin == DENSITY::Total) occ = occ_a + occ_b + occ_p;
    if (spin == DENSITY::Alpha) occ = occ_a + 0.5*occ_p;
    if (spin == DENSITY::Beta)  occ = occ_b + 0.5*occ_p;
    if (spin == DENSITY::Spin)  occ = occ_a - occ_b;

    if (std::abs(occ) < mrcpp::MachineZero) {
        rho.setZero();
        return;
    }

    // get total orbital in case there's a correlation factor
    Orbital Rphi = phi;
    if (R != nullptr) Rphi = (*R)(phi);

    FunctionTreeVector<3> sum_vec;
    if (Rphi.hasReal()) {
        FunctionTree<3> *real_2 = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*real_2, rho);
        mrcpp::multiply(prec, *real_2, occ, Rphi.real(), Rphi.real());
        sum_vec.push_back(std::make_tuple(1.0, real_2));
    }
    if (Rphi.hasImag()) {
        FunctionTree<3> *imag_2 = new FunctionTree<3>(*MRA);
        mrcpp::copy_grid(*imag_2, rho);
        mrcpp::multiply(prec, *imag_2, occ, Rphi.imag(), Rphi.imag());
        sum_vec.push_back(std::make_tuple(1.0, imag_2));
    }
    if (R != nullptr) Rphi.free();
    mrcpp::build_grid(rho, sum_vec);
    mrcpp::add(-1.0, rho, sum_vec, 0);
    mrcpp::clear(sum_vec, true);
}

void density::compute(double prec,
                      Density &rho,
                      OrbitalVector &Phi,
                      int spin,
                      NuclearCorrelationOperator *R) {
    double mult_prec = prec;            // prec for \rho_i = |\phi_i|^2
    double add_prec = prec/Phi.size();  // prec for \sum_i \rho_i

    DensityVector dens_vec;
    for (int i = 0; i < Phi.size(); i++) {
        Density *rho_i = new Density(*MRA);
        mrcpp::copy_grid(*rho_i, rho);
        density::compute(mult_prec, *rho_i, Phi[i], spin, R);
        dens_vec.push_back(std::make_tuple(1.0, rho_i));
    }

    // Adaptive prec addition if more than 5 contributions,
    // otherwise addition on union grid
    if (dens_vec.size() > 5 and add_prec > 0.0) {
        mrcpp::add(add_prec, rho, dens_vec);
    } else if (dens_vec.size() > 0) {
        mrcpp::build_grid(rho, dens_vec);
        mrcpp::add(-1.0, rho, dens_vec, 0);
    }
    mrcpp::clear(dens_vec, true);
}

} //namespace mrchem
