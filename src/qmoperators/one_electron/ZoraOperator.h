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
    ZoraOperator(QMPotential &V, double c, std::shared_ptr<mrcpp::DerivativeOperator<3>> D, bool mpi_share = false) 
        : light_speed(c)
        , derivative(D) {
            // Construct pointers
            auto vz = std::make_shared<QMPotential>(1, mpi_share);
            auto vz_inv = std::make_shared<QMPotential>(1, mpi_share);
            auto base_2cc = std::make_shared<QMPotential>(1, mpi_share);
            
            double twocc = this->two_cc();
            
            // Compute potentials that depend on base
            qmfunction::deep_copy(*vz, V);
            qmfunction::deep_copy(*base_2cc, V);
            
            auto map_vz = [twocc](double val) -> double { return twocc / (twocc - val); };
            auto map_base_over_2cc = [twocc](double val) -> double { val / twocc; };
            
            vz->real().map(map_vz);
            base_2cc->real().map(map_base_over_2cc);
            
            // Compute inverse kappa potential
            qmfunction::deep_copy(*vz_inv, *vz);
            auto map_vz_inv = [twocc](double val) -> double { return 1.0 / val; };
            vz_inv->real().map(map_vz_inv);
            
            // Assign potentials to members
            this->kappa_inv = vz_inv;
            this->base_over_2cc = base_2cc;
            
            // Invoke operator= to assign *this operator
            RankZeroOperator &O = (*this);
            O = vz;
            O.name() = "V_zora";
            
            // Compute and assign gradient of kappa
            NablaOperator nabla(this->derivative);
            RankOneOperator<3> kappa_grad = nabla(*this);
            this->kappa_grad = std::make_shared<RankOneOperator<3>>(kappa_grad);
            
            // Set operator names
            this->kappa_inv->real().setName("ZORA_inverse");
            this->base_over_2cc->real().setName("ZORA_base_over_2cc");
    };
    
public:
    double light_speed;
    std::shared_ptr<mrcpp::DerivativeOperator<3>> derivative{nullptr};
    std::shared_ptr<RankOneOperator<3>> kappa_grad;
    std::shared_ptr<QMPotential> kappa_inv{nullptr};
    std::shared_ptr<QMPotential> base_over_2cc{nullptr}; 
    
public:
    double two_cc() { return 2.0 * this->light_speed * this->light_speed; }
    RankZeroOperator kappa_pot_inv() { return RankZeroOperator(this->kappa_inv); }
    RankZeroOperator base_pot_over_2cc() { return RankZeroOperator(this->base_over_2cc); }
    RankOneOperator<3> kappa_pot_grad() { return *(this->kappa_grad); }
 
};
   
} // namespace mrchem
