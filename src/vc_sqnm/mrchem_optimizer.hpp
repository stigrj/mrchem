#include <iostream>
#include <string>

#include <MRCPP/Timer>
#include <MRCPP/Parallel>

#include "driver.h"
#include "mrchem.h"
#include "mrenv.h"
#include "version.h"

#include "chemistry/Molecule.h"
// #include "chemistry/PhysicalConstants.h"
#include "vc_sqnm/periodic_optimizer.hpp"

#include <Eigen/Dense>

using json = nlohmann::json;
using namespace mrchem;

/**
 * @brief Gets postion from a molecule in json format.
 * @param mol_inp Molecule in json format.
 * 
 * @return Positions as an Eigen::matrix of shape (3, num_atoms)
*/
Eigen::MatrixXd getPositions(const json &mol_inp) {
    Molecule mol;
    driver::init_molecule(mol_inp, mol);
    Eigen::MatrixXd pos(3, mol.getNNuclei());

    int i = 0;
    for (const auto &coord : mol_inp["coords"]) {
        pos.col(i) << coord["xyz"][0], coord["xyz"][1], coord["xyz"][2];
        i++;
    }
    return pos;
}

/**
 * @brief Sets positions in json molecule.
 * @param mol_inp Inpunt molecule in json format.
 * @param pos Eigen::MatrixXd containing positions. Dimension must be (3, num_atoms)
*/
void setPositions(json &mol_inp, const Eigen::MatrixXd &pos) {
    for (int i = 0; i < mol_inp["coords"].size(); i++) {
        Eigen::VectorXd vvv = pos.col(i);
        mol_inp["coords"][i]["xyz"] = {vvv(0), vvv(1), vvv(2)};
    }
}

/**
 * @brief Does an scf calculation of a molecule.
 * 
 * @param mol_inp: json that describes the molecule.
 * @param scf_inp: scf settings.
 * 
 * @return: json summary of scf results.
*/
json getSCFResults(const json mol_inp, const json scf_inp) {
    Molecule mol;
    driver::init_molecule(mol_inp, mol);
    auto scf_out = driver::scf::run(scf_inp, mol);
    // TODO: Ask someone if this mpi barrier is needed!
    mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk);
    return driver::print_properties(mol);
}

/**
 * @brief Extracts forces from json containing scf results.
 * @param scf_results json containing the results of an scf calculations. Must be obtained using the function getSCFResults.
 * @param forces: Provide empty matrix of dimension (3, num_atoms).
 * 
 * @return Forces acting on nuclei
*/
Eigen::MatrixXd extractForcesInPlace(const json scf_results, Eigen::MatrixXd &forces){
    for (int i = 0; i < forces.cols(); i++)
    {
        for (int j = 0; j < 3; j++) {
            // is this safe? will the forces always be stored in ["geometric_derivative"]["geom-1"]?
            forces(j, i) = scf_results["geometric_derivative"]["geom-1"]["total"][3*i + j];
        }
    }
    // forces are the negative nuclear gradient.
    return -forces;
}

/**
 * @brief Extracts energy from json containing scf results.
 * @param scf_results json containing the results of an scf calculations. Must be obtained using the function getSCFResults.
 * 
 * @return Total energy
*/
double extractEnergy(const json scf_results){
    return scf_results["scf_energy"]["E_tot"];
}

/**
 * @brief Optimizes positions of nuclei.
 * 
 * @param scf_inp: scf settings.
 * @param mol_inp: json that contains the molecule.
 * @param geopt_inp: json that contains the geometry optization settings.
 * @return: A summary of the geometry optimization trajectory.
*/
json optimize_positions(json scf_inp, json mol_inp, json geopt_inp) {

    int num_atoms = mol_inp["coords"].size();

    // define default parameters
    // The sqnm parameters are documented in the periodic_optimizer.hpp file.
    auto max_iter = geopt_inp["max_iter"];
    auto max_history_length = geopt_inp["max_history_length"];
    auto init_step_size = geopt_inp["init_step_size"];
    auto minimal_step_size = geopt_inp["minimal_step_size"];
    auto subspace_tolerance = geopt_inp["subspace_tolerance"];
    auto max_force_component = geopt_inp["max_force_component"];

    PES_optimizer::periodic_optimizer optimizer(num_atoms, init_step_size, max_history_length, minimal_step_size, subspace_tolerance);
    
    int i = 0;
    Eigen::MatrixXd forces(3, num_atoms);
    double energy;
    
    json results = getSCFResults(mol_inp, scf_inp);
    energy = extractEnergy(results);
    forces = extractForcesInPlace(results, forces);
    // write progress to stderr for earier debugging.
    std::cerr << i << " " << energy << " " << forces.cwiseAbs().maxCoeff() << "\n";

    json summary;
    summary[i] = {
        {"results", results},
        {"molecule", mol_inp}
    };

    Eigen::MatrixXd pos = getPositions(mol_inp);

    while (i < max_iter && forces.cwiseAbs().maxCoeff() > max_force_component) {
        optimizer.step(pos, energy, forces);
        setPositions(mol_inp, pos);
        if (geopt_inp["use_previous_guess"]) {
            scf_inp["initial_guess"]["type"] = "mw";
            scf_inp["initial_guess"]["file_phi_p"] = scf_inp["write_orbitals"]["file_phi_p"];
            scf_inp["initial_guess"]["file_phi_a"] = scf_inp["write_orbitals"]["file_phi_a"];
            scf_inp["initial_guess"]["file_phi_b"] = scf_inp["write_orbitals"]["file_phi_b"];
        }
        json results = getSCFResults(mol_inp, scf_inp);
        energy = extractEnergy(results);
        forces = extractForcesInPlace(results, forces);
        i++;
        summary[i] = {
            {"results", results},
            {"molecule", mol_inp}
        };
        std::cerr << i << " " << energy << " " << forces.cwiseAbs().maxCoeff() << "\n";
    }
    std::cerr << "Estimated energy of the current local minimum:                        " << optimizer.lower_bound() << "\n";
    std::cerr << "Energy difference of current energy to estimated ground state energy: " << energy - optimizer.lower_bound() << "\n";
    return summary;

}
