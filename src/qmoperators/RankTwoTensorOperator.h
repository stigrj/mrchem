#pragma once

#include "TensorOperator.h"
#include "RankOneTensorOperator.h"

namespace mrchem {

/** @class RankTwoTensorOperator
 *
 *  @brief Matrix of RankZeroTensorOperator
 *
 * This class provides a base for all matrix operators (implemented as a vector
 * of vectors), and implements some simple collective operations returning matrix
 * quantities.
 *
 */

template<int I, int J>
class RankTwoTensorOperator : public TensorOperator<I, RankOneTensorOperator<J> > {
public:
    ComplexMatrix operator()(Orbital bra, Orbital ket) {
        RankTwoTensorOperator &O = *this;
        ComplexMatrix out(I,J);
        for (int i = 0; i < I; i++) {
            out.row(i) = O[i](bra, ket);
        }
        return out;
    }
    ComplexMatrix trace(OrbitalVector &phi, NuclearCorrelationOperator *R = nullptr) {
        RankTwoTensorOperator &O = *this;
        ComplexMatrix out(I,J);
        for (int i = 0; i < I; i++) {
            out.row(i) = O[i].trace(phi, R);
        }
        return out;
    }
    ComplexMatrix trace(OrbitalVector &phi,
                        OrbitalVector &x,
                        OrbitalVector &y,
                        NuclearCorrelationOperator *R = nullptr) {
        RankTwoTensorOperator &O = *this;
        ComplexMatrix out(I,J);
        for (int i = 0; i < I; i++) {
            out.row(i) = O[i].trace(phi, x, y, R);
        }
        return out;
    }
};

} //namespace mrchem
