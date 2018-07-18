#pragma once

#include "QMOperator.h"

namespace mrchem {

/** @class RankZeroTensorOperator
 *
 *  @brief General expansion of QMOperators
 *
 * This class represents a general sum of products of fundamental QM operators,
 * which will be applied from right to left:
 *      \hat{Q} = c_1 \hat{q}_{11}\hat{q}_{12}\cdots
 *              + c_2 \hat{q}_{21}\hat{q}_{22}\cdots
 *              + \cdots
 *
 * The class has two main purposes:
 * (1)  it provides an easy, intuitive way to contruct new (more complicated) operators from
 *      simple fundamental operators by implementing the arithmetic operators (+,-,*). E.g.
 *      the x component of the angular momentum can be written explicitly from position
 *      and linear momentum as: l_x = r_y*p_z - r_z*p_y
 * (2)  it provides the public interface to all operators by implementing the operator()
 *      function, as well as expectation value and trace operation, using the (protected)
 *      apply() function of the fundamental QMOperators
 *
 * Operators can be constructed either on the fly by combining other existing operators,
 * or by defining a new derived class that contains the fundamental QMOperators (see e.g.
 * NuclearOperator which CONTAINS a NuclearPotential that is ASSIGNED to the operator itself
 * in the constructor). Note that the assignment operator is a shallow copy that does NOT
 * transfer ownership of the QMOperator pointers, so be careful when defining new operators
 * on the fly, as the underlying operators can go out of scope. This is why all implemented
 * operators contain their fundamental operators as data members.
 *
 */

// Convenience typedef
typedef std::vector<QMOperator *> QMOperatorVector;

class Nucleus;
class Nuclei;
 
class RankZeroTensorOperator {
public:
    RankZeroTensorOperator() = default;
    RankZeroTensorOperator(QMOperator &O) { *this = O; }
    RankZeroTensorOperator(const RankZeroTensorOperator &O) { *this = O; }
    virtual ~RankZeroTensorOperator() = default;

    void setup(double prec);
    void clear();

    Orbital operator()(Orbital inp);
    Orbital dagger(Orbital inp);

    OrbitalVector operator()(OrbitalVector &inp);
    OrbitalVector dagger(OrbitalVector &inp);

    ComplexDouble operator()(Orbital bra, Orbital ket, NuclearCorrelationOperator *R = nullptr);
    ComplexDouble dagger(Orbital bra, Orbital ket, NuclearCorrelationOperator *R = nullptr);

    ComplexMatrix operator()(OrbitalVector &bra, OrbitalVector &ket, NuclearCorrelationOperator *R = nullptr);
    ComplexMatrix dagger(OrbitalVector &bra, OrbitalVector &ket, NuclearCorrelationOperator *R = nullptr);

    ComplexDouble trace(OrbitalVector &Phi, NuclearCorrelationOperator *R = nullptr);
    ComplexDouble trace(OrbitalVector &Phi, OrbitalVector &X, OrbitalVector &Y, NuclearCorrelationOperator *R = nullptr);
    
    RankZeroTensorOperator& operator=(QMOperator &O);
    RankZeroTensorOperator& operator+=(QMOperator &O);
    RankZeroTensorOperator& operator-=(QMOperator &O);
    RankZeroTensorOperator& operator=(const RankZeroTensorOperator &O);
    RankZeroTensorOperator& operator+=(const RankZeroTensorOperator &O);
    RankZeroTensorOperator& operator-=(const RankZeroTensorOperator &O);

    friend RankZeroTensorOperator operator*(ComplexDouble a, RankZeroTensorOperator A);
    friend RankZeroTensorOperator operator*(RankZeroTensorOperator A, RankZeroTensorOperator B);
    friend RankZeroTensorOperator operator+(RankZeroTensorOperator A, RankZeroTensorOperator B);
    friend RankZeroTensorOperator operator-(RankZeroTensorOperator A, RankZeroTensorOperator B);

protected:
    std::vector<ComplexDouble> coef_exp;
    std::vector<QMOperatorVector> oper_exp;

    Orbital applyOperTerm(int n, Orbital inp);
    ComplexVector getCoefVector() const;
};


inline RankZeroTensorOperator operator*(ComplexDouble a, RankZeroTensorOperator A) {
    RankZeroTensorOperator out;
    for (int i = 0; i < A.oper_exp.size(); i++) {
        out.coef_exp.push_back(a*A.coef_exp[i]);
        out.oper_exp.push_back(A.oper_exp[i]);
    }
    return out;
}

inline RankZeroTensorOperator operator*(RankZeroTensorOperator A,
                                        RankZeroTensorOperator B) {
    RankZeroTensorOperator out;
    int a_terms = A.oper_exp.size();
    int b_terms = B.oper_exp.size();

    for (int i = 0; i < a_terms; i++) {
        ComplexDouble a_i = A.coef_exp[i];
        for (int k = 0; k < b_terms; k++) {
            ComplexDouble b_k = B.coef_exp[k];
            QMOperatorVector tmp;
            for (int l = 0; l < B.oper_exp[k].size(); l++) {
                QMOperator *B_kl = B.oper_exp[k][l];
                tmp.push_back(B_kl);
            }
            for (int j = 0; j < A.oper_exp[i].size(); j++) {
                QMOperator *A_ij = A.oper_exp[i][j];
                tmp.push_back(A_ij);
            }
            out.coef_exp.push_back(b_k*a_i);
            out.oper_exp.push_back(tmp);
        }
    }
    return out;
}

inline RankZeroTensorOperator operator+(RankZeroTensorOperator A,
                                        RankZeroTensorOperator B) {
    RankZeroTensorOperator out;
    for (int i = 0; i < A.oper_exp.size(); i++) {
        out.coef_exp.push_back(A.coef_exp[i]);
        out.oper_exp.push_back(A.oper_exp[i]);
    }
    for (int i = 0; i < B.oper_exp.size(); i++) {
        out.coef_exp.push_back(B.coef_exp[i]);
        out.oper_exp.push_back(B.oper_exp[i]);
    }
    return out;
}

inline RankZeroTensorOperator operator-(RankZeroTensorOperator A,
                                        RankZeroTensorOperator B) {
    RankZeroTensorOperator out;
    for (int i = 0; i < A.oper_exp.size(); i++) {
        out.coef_exp.push_back(A.coef_exp[i]);
        out.oper_exp.push_back(A.oper_exp[i]);
    }
    for (int i = 0; i < B.oper_exp.size(); i++) {
        out.coef_exp.push_back(-B.coef_exp[i]);
        out.oper_exp.push_back(B.oper_exp[i]);
    }
    return out;
}

} //namespace mrchem
