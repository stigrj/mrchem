/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2021 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#pragma once

#include "MRCPP/MWOperators"

#include "tensor/RankZeroOperator.h"
#include "qmfunctions/qmfunction_utils.h"
#include "qmoperators/QMPotential.h"
#include "qmoperators/one_electron/NablaOperator.h"

namespace mrchem {

class ZoraOperator final : public RankZeroOperator {
public:
    ZoraOperator(double c, std::shared_ptr<mrcpp::DerivativeOperator<3>> D)
            : light_speed(c)
            , derivative(D) {
    }
    
    RankZeroOperator kappaPotential() { return this->kappa; }
    RankZeroOperator kappaPotentialInverse() { return this->kappa_inv; }
    RankZeroOperator basePotentialOver2cc() { return this->base_over_2cc; }
    RankOneOperator<3> kappaPotentialGradient() { return this->kappa_grad; }
    
    std::shared_ptr<mrcpp::DerivativeOperator<3>> getDerivative() { return this->derivative; }

    double getLightSpeed() { return this->light_speed; }
    double getTwoCC() { return 2.0 * this->light_speed * this->light_speed; }
    
    void setupPotentials(QMPotential &vz, double prec) {
        double twocc = getTwoCC();
        
        // kappa potential
        auto map_k = [twocc](double val) -> double { return twocc / (twocc - val); };
        auto k = std::make_shared<QMPotential>(1, false);
        qmfunction::deep_copy(*k, vz);
        if (k->hasReal()) k->real().map(map_k);
        this->kappa = RankZeroOperator(k);
        this->kappa.setup(prec);
        
        // inverse kappa potential
        auto map_kinv = [](double val) -> double { return 1.0 / val; };
        auto k_inv = std::make_shared<QMPotential>(1, false);
        qmfunction::deep_copy(*k_inv, *k);
        if (k_inv->hasReal()) k_inv->real().map(map_kinv);
        this->kappa_inv = RankZeroOperator(k_inv);
        this->kappa_inv.setup(prec);
        
        // kappa gradient potential
        NablaOperator nabla(this->derivative);
        RankZeroOperator Kappa(this->kappa);
        this->kappa_grad = nabla(Kappa);
        this->kappa_grad.setup(prec);
        
        // base / 2cc potential
        auto map_vz_over_2cc = [twocc](double val) -> double { return val / twocc; };
        auto vz_over_2cc = std::make_shared<QMPotential>(1, false);
        qmfunction::deep_copy(*vz_over_2cc, vz);
        if (vz_over_2cc->hasReal()) vz_over_2cc->real().map(map_vz_over_2cc);
        this->base_over_2cc = RankZeroOperator(vz_over_2cc);
        this->base_over_2cc.setup(prec);
    }

    void clearPotentials() {
        this->kappa.clear();
        this->kappa_inv.clear();
        this->kappa_grad.clear();
        this->base_over_2cc.clear();
        this->kappa = RankZeroOperator();
        this->kappa_inv = RankZeroOperator();
        this->kappa_grad = RankOneOperator<3>();
        this->base_over_2cc = RankZeroOperator();
    }

private:
    double light_speed;
    std::shared_ptr<mrcpp::DerivativeOperator<3>> derivative{nullptr};
    
    RankZeroOperator kappa;
    RankZeroOperator kappa_inv;
    RankZeroOperator base_over_2cc;
    RankOneOperator<3> kappa_grad;

};
   
} // namespace mrchem
