/*
 * MRChem, a numerical real-space code for molecular electronic structure
 * calculations within the self-consistent field (SCF) approximations of quantum
 * chemistry (Hartree-Fock and Density Functional Theory).
 * Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
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

#include "MRCPP/MWOperators"
#include "MRCPP/Printer"
#include "MRCPP/Timer"
#include <MRCPP/trees/FunctionNode.h>

#include "Functional.h"
#include "MRDFT.h"
#include "utils/Bank.h"
#include "xc_utils.h"

namespace mrdft {

/** @brief Compute XC potentials from densities
 *
 * This routine computes the XC energy density and potentials on the
 * union grid of all the density input functions. Functional evaluation
 * and subsequent contraction is done node by node, to avoid explicit
 * construction of the huge number of intermediate functions.
 *
 * Ordering without spin:
 * inp_vec[0] = rho_0 (unperturbed)
 * inp_vec[1] = rho_1 (first order perturbed)
 * ...
 * out_vec[0] = f_xc (XC energy density)
 * out_vec[1] = v_xc (XC potential)
 *
 * Ordering with spin:
 * inp_vec[0] = alpha_0 (unperturbed)
 * inp_vec[1] = beta_0
 * inp_vec[2] = alpha_1 (first order perturbed)
 * inp_vec[3] = beta_1
 * ...
 * out_vec[0] = f_xc (XC energy density)
 * out_vec[1] = v_xc_a (XC alpha potential)
 * out_vec[2] = v_xc_b (XC beta potential)
 */
mrcpp::FunctionTreeVector<3> MRDFT::evaluate(mrcpp::FunctionTreeVector<3> &inp) {
    mrcpp::Timer t_tot, t_pre;
    grid().unify(inp);
    functional().preprocess(inp);
    mrcpp::FunctionTreeVector<3> xcInpVec = functional().setupXCInput();
    mrcpp::FunctionTreeVector<3> ctrInpVec = functional().setupCtrInput();

    int inpNodes = 0;
    int inpSize = 0;
    for (auto i = 0; i < xcInpVec.size(); i++) {
        auto &f_i = mrcpp::get_func(xcInpVec, i);
        inpNodes += f_i.getNNodes();
        inpSize += f_i.getSizeNodes();
    }
    mrcpp::print::tree(3, "Preprocess input", inpNodes, inpSize, t_pre.elapsed());

    mrcpp::Timer t_eval;
    int nCoefs = mrcpp::get_func(inp, 0).getEndFuncNode(0).getNCoefs();
    int nOutCtr = functional().getCtrOutputLength();
    int nFcs = functional().getXCOutputLength();

    // divide nNodes into parts assigned to each MPI rank
    int nNodes = grid().size();
    int n_start = (mrchem::mpi::orb_rank * nNodes) / mrchem::mpi::orb_size;
    int n_end = ((mrchem::mpi::orb_rank + 1) * nNodes) / mrchem::mpi::orb_size;
    std::vector<Eigen::MatrixXd> ctrOutDataVec(n_end - n_start);

    mrcpp::FunctionTreeVector<3> ctrOutVec;
    if (mrchem::mpi::orb_size == 1) ctrOutVec = grid().generate(nOutCtr);

#pragma omp parallel
    {
#pragma omp for schedule(guided)
        for (int n = n_start; n < n_end; n++) {
            auto xcInpNodes = xc_utils::fetch_nodes(n, xcInpVec);
            auto xcInpData = xc_utils::compress_nodes(xcInpNodes);

            auto xcOutData = functional().evaluate(xcInpData);
            auto ctrInpNodes = xc_utils::fetch_nodes(n, ctrInpVec);
            auto ctrInpData = xc_utils::compress_nodes(ctrInpNodes);
            auto ctrOutData = functional().contract(xcOutData, ctrInpData);

            if (mrchem::mpi::orb_size > 1) {
                // store the results temporarily
                ctrOutDataVec[n - n_start] = std::move(ctrOutData);
            } else {
                // postprocess the results
                auto ctrOutNodes = xc_utils::fetch_nodes(n, ctrOutVec);
                xc_utils::expand_nodes(ctrOutNodes, ctrOutData);
            }
        }
    }

    // Input data is cleared before constructing the full output
    mrcpp::clear(xcInpVec, false);
    mrcpp::clear(ctrInpVec, false);

    if (mrchem::mpi::orb_size > 1) {
        // each MPI process has only a part of the results

        ctrOutVec = grid().generate(nOutCtr);
        mrchem::BankAccount ctrOutBank; // to put the ctrOutDataVec;

        // note that mpi cannot run in multiple omp threads
        int size = nOutCtr * nCoefs;
        for (int n = n_start; n < n_end; n++) ctrOutBank.put_data(n, size, ctrOutDataVec[n - n_start].data());
        // fetch all nodes from bank and postprocess
        for (int n = 0; n < nNodes; n++) {
            Eigen::MatrixXd ctrOutData(nOutCtr, nCoefs);
            ctrOutBank.get_data(n, size, ctrOutData.data());
            auto ctrOutNodes = xc_utils::fetch_nodes(n, ctrOutVec);
            xc_utils::expand_nodes(ctrOutNodes, ctrOutData);
        }
    }

    // Reconstruct raw xcfun output functions
    /*
    for (auto i = 0; i < xcOutVec.size(); i++) {
        auto &f_i = mrcpp::get_func(xcOutVec, i);
        f_i.mwTransform(mrcpp::BottomUp);
        f_i.calcSquareNorm();
    }
    mrcpp::clear(xcOutVec, true);
    */

    // Reconstruct contracted output functions
    int ctrNodes = 0;
    int ctrSize = 0;
    for (auto i = 0; i < ctrOutVec.size(); i++) {
        auto &f_i = mrcpp::get_func(ctrOutVec, i);
        f_i.mwTransform(mrcpp::BottomUp);
        f_i.calcSquareNorm();
        // TODO? insert a crop
        ctrNodes += f_i.getNNodes();
        ctrSize += f_i.getSizeNodes();
    }
    mrcpp::print::tree(3, "Evaluate functional", ctrNodes, ctrSize, t_eval.elapsed());

    mrcpp::Timer t_post;
    auto potOutVec = functional().postprocess(ctrOutVec);
    mrcpp::clear(ctrOutVec, true);
    functional().clear();

    int outNodes = 0;
    int outSize = 0;
    for (auto i = 0; i < potOutVec.size(); i++) {
        auto &f_i = mrcpp::get_func(potOutVec, i);
        f_i.mwTransform(mrcpp::BottomUp);
        f_i.calcSquareNorm();
        // TODO? insert a crop
        outNodes += f_i.getNNodes();
        outSize += f_i.getSizeNodes();
    }
    mrcpp::print::tree(3, "Postprocess potential", outNodes, outSize, t_post.elapsed());

    return potOutVec;
}

} // namespace mrdft
