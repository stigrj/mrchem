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

#include "catch.hpp"

#include <MRCPP/MWOperators>
#include <MRCPP/Printer>

#include "mrchem.h"
#include "parallel.h"

#include "analyticfunctions/HarmonicOscillatorFunction.h"
#include "qmfunctions/Orbital.h"
#include "qmfunctions/orbital_utils.h"
#include "qmoperators/one_electron/IdentityOperator.h"
#include "qmoperators/one_electron/MomentumOperator.h"
#include "qmoperators/one_electron/NablaOperator.h"
#include "qmoperators/one_electron/PositionOperator.h"
#include "qmoperators/one_electron/SpinOperator.h"
#include "tensor/RankTwoOperator.h"
#include "tensor/tensor_utils.h"

using namespace mrchem;

namespace operator_composition {

TEST_CASE("Operator composition", "[operator_composition]") {
    const double prec = 1.0e-3;
    const double thrs = prec * prec;

    Orbital phi_0(SPIN::Alpha), phi_1(SPIN::Alpha), phi_2(SPIN::Alpha);
    {
        int nu[3] = {0, 0, 0};
        HarmonicOscillatorFunction f(nu);
        qmfunction::project(phi_0, f, NUMBER::Real, prec);
    }
    {
        int nu[3] = {1, 0, 0};
        HarmonicOscillatorFunction f(nu);
        qmfunction::project(phi_1, f, NUMBER::Real, prec);
    }
    {
        int nu[3] = {2, 0, 0};
        HarmonicOscillatorFunction f(nu);
        qmfunction::project(phi_2, f, NUMBER::Real, prec);
    }

    SECTION("identity operator") {
        IdentityOperator I;
        SECTION("product I*I") {
            RankZeroOperator II = I * I;
            REQUIRE(II.size() == 1);
            REQUIRE(II.size(0) == 2);

            II.setup(prec);
            const double ref = phi_0.squaredNorm();
            const ComplexDouble val = II(phi_0, phi_0);
            REQUIRE(val.real() == Approx(ref).epsilon(thrs));
            REQUIRE(val.imag() == Approx(0.0).margin(thrs));
            II.clear();
        }
        SECTION("composition I(I)") {
            RankZeroOperator II = I(I);
            REQUIRE(II.size() == 1);
            REQUIRE(II.size(0) == 1);

            II.setup(prec);
            const double ref = phi_0.squaredNorm();
            const ComplexDouble val = II(phi_0, phi_0);
            REQUIRE(val.real() == Approx(ref).epsilon(thrs));
            REQUIRE(val.imag() == Approx(0.0).margin(thrs));
            II.clear();
        }
        SECTION("product I*V") {
            PositionOperator V;
            RankZeroOperator IV = I * V[0];
            REQUIRE(IV.size() == 1);
            REQUIRE(IV.size(0) == 2);

            IV.setup(prec);
            const ComplexDouble ref = V[0](phi_1, phi_0);
            const ComplexDouble val = IV(phi_1, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).epsilon(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).margin(thrs));
            IV.clear();
        }
        SECTION("composition I(V)") {
            PositionOperator V;
            RankZeroOperator IV = I(V[0]);
            REQUIRE(IV.size() == 1);
            REQUIRE(IV.size(0) == 1);

            IV.setup(prec);
            const ComplexDouble ref = V[0](phi_1, phi_0);
            const ComplexDouble val = IV(phi_1, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).epsilon(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).margin(thrs));
            IV.clear();
        }
        SECTION("product I*D") {
            MomentumOperator D(std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5));
            RankZeroOperator ID = I * D[0];
            REQUIRE(ID.size() == 1);
            REQUIRE(ID.size(0) == 2);

            ID.setup(prec);
            const ComplexDouble ref = D[0](phi_0, phi_1);
            const ComplexDouble val = ID(phi_0, phi_1);
            REQUIRE(val.real() == Approx(ref.real()).margin(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).epsilon(thrs));
            ID.clear();
        }
        SECTION("composition I(D)") {
            MomentumOperator D(std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5));
            RankZeroOperator ID = I(D[0]);
            REQUIRE(ID.size() == 1);
            REQUIRE(ID.size(0) == 1);

            ID.setup(prec);
            const ComplexDouble ref = D[0](phi_0, phi_1);
            const ComplexDouble val = ID(phi_0, phi_1);
            REQUIRE(val.real() == Approx(ref.real()).margin(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).epsilon(thrs));
            ID.clear();
        }
        SECTION("product I*S") {
            SpinOperator S;
            RankZeroOperator IS = I * S[1];
            REQUIRE(IS.size() == 1);
            REQUIRE(IS.size(0) == 2);

            IS.setup(prec);
            const ComplexDouble ref = S[1](phi_0, phi_0);
            const ComplexDouble val = IS(phi_0, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).margin(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).epsilon(thrs));
            IS.clear();
        }
        SECTION("composition I(S)") {
            SpinOperator S;
            RankZeroOperator IS = I(S[1]);
            REQUIRE(IS.size() == 1);
            REQUIRE(IS.size(0) == 1);

            IS.setup(prec);
            const ComplexDouble ref = S[1](phi_0, phi_0);
            const ComplexDouble val = IS(phi_0, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).margin(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).epsilon(thrs));
            IS.clear();
        }
    }
    SECTION("potential operator") {
        PositionOperator V;
        SECTION("product V*I") {
            IdentityOperator I;
            RankZeroOperator VI = V[0] * I;
            REQUIRE(VI.size() == 1);
            REQUIRE(VI.size(0) == 2);

            VI.setup(prec);
            const ComplexDouble ref = {std::sqrt(2.0) / 2.0, 0.0};
            const ComplexDouble val = VI(phi_1, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).epsilon(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).margin(thrs));
            VI.clear();
        }
        SECTION("composition V(I)") {
            IdentityOperator I;
            RankZeroOperator VI = V[0](I);
            REQUIRE(VI.size() == 1);
            REQUIRE(VI.size(0) == 1);

            VI.setup(prec);
            const ComplexDouble ref = {std::sqrt(2.0) / 2.0, 0.0};
            const ComplexDouble val = VI(phi_1, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).epsilon(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).margin(thrs));
            VI.clear();
        }
        SECTION("product V*V") {
            RankZeroOperator VV = V[0] * V[0];
            REQUIRE(VV.size() == 1);
            REQUIRE(VV.size(0) == 2);

            VV.setup(prec);
            const ComplexDouble ref = {std::sqrt(2.0) / 2.0, 0.0};
            const ComplexDouble val = VV(phi_2, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).epsilon(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).margin(thrs));
            VV.clear();
        }
        SECTION("composition V(V)") {
            RankZeroOperator VV = V[0](V[0]);
            REQUIRE(VV.size() == 1);
            REQUIRE(VV.size(0) == 1);

            VV.setup(prec);
            const ComplexDouble ref = {std::sqrt(2.0) / 2.0, 0.0};
            const ComplexDouble val = VV(phi_2, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).epsilon(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).margin(thrs));
            VV.clear();
        }
        SECTION("product V*D") {
            MomentumOperator D(std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5));
            RankZeroOperator VD = V[0] * D[0];
            REQUIRE(VD.size() == 1);
            REQUIRE(VD.size(0) == 2);

            VD.setup(prec);
            const ComplexDouble ref = {0.0, 0.5};
            const ComplexDouble val = VD(phi_0, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).margin(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).epsilon(thrs));
            VD.clear();
        }
        SECTION("composition V(D)") {
            MomentumOperator D(std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5));
            RankZeroOperator VD = V[0](D[0]);
            REQUIRE(VD.size() == 1);
            REQUIRE(VD.size(0) == 2);

            VD.setup(prec);
            const ComplexDouble ref = {0.0, 0.5};
            const ComplexDouble val = VD(phi_0, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).margin(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).epsilon(thrs));
            VD.clear();
        }
        SECTION("product V*S") {
            SpinOperator S;
            RankZeroOperator VS = V[0] * S[2];
            REQUIRE(VS.size() == 1);
            REQUIRE(VS.size(0) == 2);

            VS.setup(prec);
            const ComplexDouble ref = {std::sqrt(2.0) / 4.0, 0.0};
            const ComplexDouble val = VS(phi_1, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).epsilon(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).margin(thrs));
            VS.clear();
        }
        SECTION("composition V(S)") {
            SpinOperator S;
            RankZeroOperator VS = V[0](S[2]);
            REQUIRE(VS.size() == 1);
            REQUIRE(VS.size(0) == 2);

            VS.setup(prec);
            const ComplexDouble ref = {std::sqrt(2.0) / 4.0, 0.0};
            const ComplexDouble val = VS(phi_1, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).epsilon(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).margin(thrs));
            VS.clear();
        }
    }
    SECTION("derivative operator") {
        MomentumOperator D(std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5));
        SECTION("product D*I") {
            IdentityOperator I;
            RankZeroOperator DI = D[0] * I;
            REQUIRE(DI.size() == 1);
            REQUIRE(DI.size(0) == 2);

            DI.setup(prec);
            const ComplexDouble ref = {0.0, std::sqrt(2.0) / 2.0};
            const ComplexDouble val = DI(phi_1, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).margin(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).epsilon(thrs));
            DI.clear();
        }
        SECTION("composition D(I)") {
            IdentityOperator I;
            RankZeroOperator DI = D[0](I);
            REQUIRE(DI.size() == 1);
            REQUIRE(DI.size(0) == 1);

            DI.setup(prec);
            const ComplexDouble ref = {0.0, std::sqrt(2.0) / 2.0};
            const ComplexDouble val = DI(phi_1, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).margin(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).epsilon(thrs));
            DI.clear();
        }
        SECTION("product D*V") {
            PositionOperator V;
            RankZeroOperator DV = D[0] * V[0];
            REQUIRE(DV.size() == 1);
            REQUIRE(DV.size(0) == 2);

            DV.setup(prec);
            const ComplexDouble ref = {0.0, -0.5};
            const ComplexDouble val = DV(phi_0, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).margin(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).epsilon(thrs));
            DV.clear();
        }
        SECTION("composition D(V)") {
            PositionOperator V;
            RankZeroOperator DV = D[0](V[0]);
            REQUIRE(DV.size() == 1);
            REQUIRE(DV.size(0) == 1);

            DV.setup(prec);
            const ComplexDouble ref = {0.0, -1.0};
            const ComplexDouble val = DV(phi_0, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).margin(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).epsilon(thrs));
            DV.clear();
        }
        SECTION("product D*D") {
            RankZeroOperator DD = D[0] * D[0];
            REQUIRE(DD.size() == 1);
            REQUIRE(DD.size(0) == 2);

            DD.setup(prec);
            const ComplexDouble ref = {0.5, 0.0};
            const ComplexDouble val = DD(phi_0, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).epsilon(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).margin(thrs));
            DD.clear();
        }
        SECTION("composition D(D)") {
            RankZeroOperator DD = D[0](D[0]);
            REQUIRE(DD.size() == 1);
            REQUIRE(DD.size(0) == 2);

            DD.setup(prec);
            const ComplexDouble ref = {0.5, 0.0};
            const ComplexDouble val = DD(phi_0, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).epsilon(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).margin(thrs));
            DD.clear();
        }
        SECTION("product D*S") {
            SpinOperator S;
            RankZeroOperator DS = D[0] * S[2];
            REQUIRE(DS.size() == 1);
            REQUIRE(DS.size(0) == 2);

            DS.setup(prec);
            const ComplexDouble ref = {0.0, std::sqrt(2.0) / 4.0};
            const ComplexDouble val = DS(phi_1, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).margin(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).epsilon(thrs));
            DS.clear();
        }
        SECTION("composition D(S)") {
            SpinOperator S;
            RankZeroOperator DS = D[0](S[2]);
            REQUIRE(DS.size() == 1);
            REQUIRE(DS.size(0) == 2);

            DS.setup(prec);
            const ComplexDouble ref = {0.0, std::sqrt(2.0) / 4.0};
            const ComplexDouble val = DS(phi_1, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).margin(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).epsilon(thrs));
            DS.clear();
        }
    }
    SECTION("spin operator") {
        SpinOperator S;
        SECTION("product S*I") {
            IdentityOperator I;
            RankZeroOperator SI = S[1] * I;
            REQUIRE(SI.size() == 1);
            REQUIRE(SI.size(0) == 2);

            SI.setup(prec);
            const ComplexDouble ref = {0.0, 0.5};
            const ComplexDouble val = SI(phi_0, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).margin(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).epsilon(thrs));
            SI.clear();
        }
        SECTION("composition S(I)") {
            IdentityOperator I;
            RankZeroOperator SI = S[1](I);
            REQUIRE(SI.size() == 1);
            REQUIRE(SI.size(0) == 1);

            SI.setup(prec);
            const ComplexDouble ref = {0.0, 0.5};
            const ComplexDouble val = SI(phi_0, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).margin(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).epsilon(thrs));
            SI.clear();
        }
        SECTION("product S*V") {
            PositionOperator V;
            RankZeroOperator SV = S[2] * V[0];
            REQUIRE(SV.size() == 1);
            REQUIRE(SV.size(0) == 2);

            SV.setup(prec);
            const ComplexDouble ref = {std::sqrt(2.0) / 4.0, 0.0};
            const ComplexDouble val = SV(phi_1, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).epsilon(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).margin(thrs));
            SV.clear();
        }
        SECTION("product S(V)") {
            PositionOperator V;
            RankZeroOperator SV = S[2](V[0]);
            REQUIRE(SV.size() == 1);
            REQUIRE(SV.size(0) == 2);

            SV.setup(prec);
            const ComplexDouble ref = {std::sqrt(2.0) / 4.0, 0.0};
            const ComplexDouble val = SV(phi_1, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).epsilon(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).margin(thrs));
            SV.clear();
        }
        SECTION("product S*D") {
            MomentumOperator D(std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5));
            RankZeroOperator SD = S[0] * D[0];
            REQUIRE(SD.size() == 1);
            REQUIRE(SD.size(0) == 2);

            SD.setup(prec);
            const ComplexDouble ref = {0.0, -std::sqrt(2.0) / 4.0};
            const ComplexDouble val = SD(phi_0, phi_1);
            REQUIRE(val.real() == Approx(ref.real()).margin(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).epsilon(thrs));
            SD.clear();
        }
        SECTION("composition S(D)") {
            MomentumOperator D(std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5));
            RankZeroOperator SD = S[0](D[0]);
            REQUIRE(SD.size() == 1);
            REQUIRE(SD.size(0) == 2);

            SD.setup(prec);
            const ComplexDouble ref = {0.0, -std::sqrt(2.0) / 4.0};
            const ComplexDouble val = SD(phi_0, phi_1);
            REQUIRE(val.real() == Approx(ref.real()).margin(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).epsilon(thrs));
            SD.clear();
        }
        SECTION("product S*S") {
            SpinOperator S;
            RankZeroOperator SS = S[1] * S[1];
            REQUIRE(SS.size() == 1);
            REQUIRE(SS.size(0) == 2);

            SS.setup(prec);
            const ComplexDouble ref = {1.0 / 4.0, 0.0};
            const ComplexDouble val = SS(phi_0, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).epsilon(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).margin(thrs));
            SS.clear();
        }
        SECTION("composition S(S)") {
            SpinOperator S;
            RankZeroOperator SS = S[1](S[1]);
            REQUIRE(SS.size() == 1);
            REQUIRE(SS.size(0) == 2);

            SS.setup(prec);
            const ComplexDouble ref = {1.0 / 4.0, 0.0};
            const ComplexDouble val = SS(phi_0, phi_0);
            REQUIRE(val.real() == Approx(ref.real()).epsilon(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).margin(thrs));
            SS.clear();
        }
    }
    SECTION("vector operators") {
        NablaOperator D(std::make_shared<mrcpp::ABGVOperator<3>>(*MRA, 0.5, 0.5));
        PositionOperator V;
        SECTION("gradient Del(V_scalar)") {
            RankOneOperator<3> gradV = D(V[0]);
            REQUIRE(gradV[0].size() == 1);
            REQUIRE(gradV[1].size() == 1);
            REQUIRE(gradV[2].size() == 1);
            REQUIRE(gradV[0].size(0) == 1);
            REQUIRE(gradV[1].size(0) == 1);
            REQUIRE(gradV[2].size(0) == 1);

            gradV.setup(prec);
            OrbitalVector dPhi_0 = gradV(phi_0);
            DoubleVector norms = orbital::get_norms(dPhi_0);
            REQUIRE(dPhi_0.size() == 3);
            REQUIRE(norms[0] == Approx(1.0).epsilon(thrs));
            REQUIRE(norms[1] == Approx(0.0).margin(thrs));
            REQUIRE(norms[2] == Approx(0.0).margin(thrs));
            gradV.clear();
        }
        SECTION("divergence De . V_vector") {
            RankZeroOperator divV = tensor::dot(D, V);
            REQUIRE(divV.size() == 3);
            REQUIRE(divV.size(0) == 1);
            REQUIRE(divV.size(1) == 1);
            REQUIRE(divV.size(2) == 1);

            divV.setup(prec);
            Orbital psi_0 = divV(phi_0);
            const ComplexDouble ref = 3.0 * phi_0.integrate();
            const ComplexDouble val = psi_0.integrate();
            REQUIRE(val.real() == Approx(ref.real()).epsilon(thrs));
            REQUIRE(val.imag() == Approx(ref.imag()).margin(thrs));
            divV.clear();
        }
        SECTION("curl Del x V_vector") {
            RankOneOperator<3> curlV = tensor::cross(D, V);
            curlV.setup(prec);
            for (int i = 0; i < 3; i++) {
                RankZeroOperator curlV_i = curlV[i];
                REQUIRE(curlV[i].size() == 2);
                REQUIRE(curlV[i].size(0) == 1);
                REQUIRE(curlV[i].size(1) == 1);
                Orbital psi_0 = curlV[i](phi_0);
                REQUIRE(psi_0.norm() == Approx(0.0).margin(thrs));
            }
            curlV.clear();
        }
        SECTION("jacobian Del (x) V_vector") {
            RankTwoOperator<3, 3> jacV = tensor::outer(D, V);
            jacV.setup(prec);
            for (int i = 0; i < 3; i++) {
                for (int j = 0; j < 3; j++) {
                    RankZeroOperator jacV_ij = jacV[i][j];
                    REQUIRE(jacV_ij.size() == 1);
                    REQUIRE(jacV_ij.size(0) == 1);
                    Orbital psi_0 = jacV_ij(phi_0);
                    if (i == j) {
                        REQUIRE(psi_0.norm() == Approx(1.0).epsilon(thrs));
                    } else {
                        REQUIRE(psi_0.norm() == Approx(0.0).margin(thrs));
                    }
                }
            }
            jacV.clear();
        }
    }
}

} // namespace operator_composition
