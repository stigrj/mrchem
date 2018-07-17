#pragma once

#include "RankZeroTensorOperator.h"

/** @class FockOperator
 *
 * @brief Operator containing the standard SCF operators
 *
 * This is a simple collection of operators used in ground state SCF calculations.
 * The operator is separated into kinetic and potential parts, since the MW way of
 * solving the SCF equations is to invert the kinetic part, and apply the potential
 * part as usual.
 */

namespace mrchem {

class SCFEnergy;
class KineticOperator;
class NuclearOperator;
class CoulombOperator;
class ExchangeOperator;
class XCOperator;
class ElectricFieldOperator;

class FockOperator final : public RankZeroTensorOperator {
public:
    FockOperator(KineticOperator        *t = 0,
                 NuclearOperator        *v = 0,
                 NuclearOperator        *u = 0,
                 CoulombOperator        *j = 0,
                 ExchangeOperator       *k = 0,
                 XCOperator             *xc = 0,
                 ElectricFieldOperator  *ext = 0);

    RankZeroTensorOperator& kinetic()       { return this->T; }
    RankZeroTensorOperator& potential()     { return this->V; }
    RankZeroTensorOperator& perturbation()  { return this->H_1; }

    KineticOperator        *getKineticOperator()     { return this->kin;  }
    NuclearOperator        *getFullNuclearOperator() { return this->nuc_full;  }
    NuclearOperator        *getNuclearOperator()     { return this->nuc_R;  }
    CoulombOperator        *getCoulombOperator()     { return this->coul; }
    ExchangeOperator       *getExchangeOperator()    { return this->ex;   }
    XCOperator             *getXCOperator()          { return this->xc;   }
    ElectricFieldOperator  *getExtOperator()         { return this->ext;  }
    
    void setKineticOperator    (KineticOperator         *t) { this->kin = t;      }
    void setFullNuclearOperator(NuclearOperator         *v) { this->nuc_full = v; }
    void setNuclearOperator    (NuclearOperator         *u) { this->nuc_R = u;    }
    void setCoulombOperator    (CoulombOperator         *j) { this->coul = j;     }
    void setExchangeOperator   (ExchangeOperator        *k) { this->ex = k;       }
    void setXCOperator         (XCOperator             *xc) { this->xc = xc;      }
    void setExtOperator        (ElectricFieldOperator *ext) { this->ext = ext;    }
    
    void rotate(const ComplexMatrix &U);

    void build();
    void setup(double prec);
    void clear();

    ComplexMatrix operator()(OrbitalVector &bra, OrbitalVector &ket, NuclearCorrelationOperator *R = nullptr);
    ComplexMatrix dagger(OrbitalVector &bra, OrbitalVector &ket, NuclearCorrelationOperator *R = nullptr);
    
    SCFEnergy trace(OrbitalVector &Phi, const ComplexMatrix &F, NuclearCorrelationOperator *R = nullptr);

    using RankZeroTensorOperator::operator();
    using RankZeroTensorOperator::dagger;

protected:
    RankZeroTensorOperator T;     ///< Total kinetic energy operator
    RankZeroTensorOperator V;     ///< Total potential energy operator
    RankZeroTensorOperator H_1;   ///< Perturbation operators

    KineticOperator        *kin;
    NuclearOperator        *nuc_full;
    NuclearOperator        *nuc_R;
    CoulombOperator        *coul;
    ExchangeOperator       *ex;
    XCOperator             *xc;
    ElectricFieldOperator  *ext;  ///< Total external potential
};

} //namespace mrchem
