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
    ZoraOperator(double c, std::shared_ptr<mrcpp::DerivativeOperator<3>> D, bool mpi_share = false)
        : light_speed(c)
        , derivative(D)
        , shared(mpi_share)
        {}
    
    double two_cc() { return 2.0 * this->light_speed * this->light_speed; }
    
    RankZeroOperator kappaPotential() { return RankZeroOperator(this->kappa); }
    RankZeroOperator kappaPotentialInverse() { return RankZeroOperator(this->kappa_inv); }
    RankZeroOperator basePotentialOver2cc() { return RankZeroOperator(this->base_over_2cc); }
    RankOneOperator<3> kappaPotentialGradient() { return *(this->kappa_grad); }
    
    std::shared_ptr<mrcpp::DerivativeOperator<3>> getDerivative() { return this->derivative; }
    double getLightSpeed() { return this->light_speed; }
    
    void updatePotentials(std::shared_ptr<RankZeroOperator> base) {
        double twocc = this->two_cc();
        auto &vz = static_cast<QMPotential &>(base->getRaw(0, 0));
        
        // kappa potential
        auto k = std::make_shared<QMPotential>(1, this->shared);
        qmfunction::deep_copy(*k, vz);
        auto map_k = [twocc](double val) -> double { return twocc / (twocc - val); };
        k->real().map(map_k);
        
        // inverse kappa potential
        auto k_inv = std::make_shared<QMPotential>(1, this->shared);
        qmfunction::deep_copy(*k_inv, *k);
        auto map_kinv = [](double val) -> double { return 1.0 / val; };
        k_inv->real().map(map_kinv);
        
        // kappa gradient potential
        NablaOperator nabla(this->derivative);
        RankOneOperator<3> k_grad = nabla(RankZeroOperator(k));
        
        // base / 2cc potential
        auto vz_over_2cc = std::make_shared<QMPotential>(1, this->shared);
        qmfunction::deep_copy(*vz_over_2cc, vz);
        auto map_vz_over_2cc = [twocc](double val) -> double { return val / twocc; };
        vz_over_2cc->real().map(map_vz_over_2cc);
        
        // Set class members
        this->kappa = k;
        this->kappa_inv = k_inv;
        this->kappa_grad = std::make_shared<RankOneOperator<3>>(k_grad);
        this->base_over_2cc = vz_over_2cc;
        }
        
    void setBasePotential(int key) {
            // Set the base potential enum from input integer
            switch (key) {
                case 0:
                    this->base_potential = NUCLEAR;
                    this->base_potential_name = "V_n";
                    break;
                case 1:
                    this->base_potential = NUCLEAR_COULOMB;
                    this->base_potential_name = "V_n + J";
                    break;
                case 2:
                    this->base_potential = NUCLEAR_COULOMB_XC;
                    this->base_potential_name = "V_n + J + V_xc";
                    break;
            }
        }
        
public:
    enum BasePotential { 
        NUCLEAR = 0, 
        NUCLEAR_COULOMB, 
        NUCLEAR_COULOMB_XC 
        };
    BasePotential base_potential;
    std::string base_potential_name;
    
private:
    double light_speed;
    bool shared;
    
private:
    std::shared_ptr<mrcpp::DerivativeOperator<3>> derivative;
    std::shared_ptr<QMPotential> kappa{nullptr};
    std::shared_ptr<QMPotential> kappa_inv{nullptr};
    std::shared_ptr<QMPotential> base_over_2cc{nullptr};
    std::shared_ptr<RankOneOperator<3>> kappa_grad{nullptr};
};
   
} // namespace mrchem
