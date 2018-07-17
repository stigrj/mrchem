#pragma once

#include "QMOperator.h"
#include "RankZeroTensorOperator.h"

namespace mrchem {

class QMIdentity final : public QMOperator {
public:
    QMIdentity() : QMOperator() { }

protected:
    void setup(double prec) { setApplyPrec(prec); }
    void clear() { clearApplyPrec(); }

    Orbital apply(Orbital inp);
    Orbital dagger(Orbital inp);
};

class IdentityOperator final : public RankZeroTensorOperator {
public:
    IdentityOperator() {
        RankZeroTensorOperator &h = (*this);
        h = I;
    }

    ComplexDouble operator()(Orbital bra, Orbital ket, NuclearCorrelationOperator *R = nullptr);
    ComplexDouble dagger(Orbital bra, Orbital ket, NuclearCorrelationOperator *R = nullptr);

    ComplexMatrix operator()(OrbitalVector &bra, OrbitalVector &ket, NuclearCorrelationOperator *R = nullptr);
    ComplexMatrix dagger(OrbitalVector &bra, OrbitalVector &ket, NuclearCorrelationOperator *R = nullptr);

    // Necessary in order to pick up base class definitions for overloaded functions
    using RankZeroTensorOperator::operator();
    using RankZeroTensorOperator::dagger;

protected:
    QMIdentity I;
};

} //namespace mrchem
