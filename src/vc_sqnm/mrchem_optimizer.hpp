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

Eigen::MatrixXd getPositions(json mol_inp) {
    Molecule mol;
    driver::init_molecule(mol_inp, mol);
    Eigen::MatrixXd pos(3, mol.getNNuclei());

    int i = 0;
    for (const auto &coord : mol_inp["coords"]) {
        double xyz[3];
        xyz[0] = coord["xyz"][0];
        xyz[1] = coord["xyz"][1];
        xyz[2] = coord["xyz"][2];

        pos.col(i) = Eigen::Vector3d(xyz);
        i++;
    }
    return pos;
}

void setPositions(json &mol_inp, Eigen::MatrixXd &pos) {
    for (int i = 0; i < mol_inp["coords"].size(); i++) {
        Eigen::VectorXd vvv = pos.col(i);
        mol_inp["coords"][i]["xyz"] = {vvv(0), vvv(1), vvv(2)};
    }
}

void energyAndForces(json mol_inp, json scf_inp, double &energy, Eigen::MatrixXd &forces) {
    Molecule mol;
    driver::init_molecule(mol_inp, mol);
    auto scf_out = driver::scf::run(scf_inp, mol);
    mrcpp::mpi::barrier(mrcpp::mpi::comm_wrk);
    json results = driver::print_properties(mol);
    energy = results["scf_energy"]["E_tot"];
    for (int i = 0; i < mol.getNNuclei(); i++)
    {
        for (int j = 0; j < 3; j++) {
            forces(j, i) = results["geometric_derivative"]["geom-1"]["total"][3*i + j];
        }
    }
    forces = -forces;
}

json optimize_positions(json scf_inp, json mol_inp, json geopt_inp) {

    int num_atoms = mol_inp["coords"].size();

    // define default parameters
    double init_step_size = -.1;
    int max_history_length = 10;
    double minimal_step_size = 0.01;
    double subspace_tolerance = 1e-3;
    int max_iter = 100;
    double max_force_component = 1e-4;

    // overwrite parameters with user defined parameters:
    if (geopt_inp.contains("init_step_size")) init_step_size = geopt_inp["init_step_size"];
    if (geopt_inp.contains("max_history_length")) max_history_length = geopt_inp["max_history_length"];
    if (geopt_inp.contains("minimal_step_size")) minimal_step_size = geopt_inp["minimal_step_size"];
    if (geopt_inp.contains("subspace_tolerance")) subspace_tolerance = geopt_inp["subspace_tolerance"];
    if (geopt_inp.contains("max_iter")) max_iter = geopt_inp["max_iter"];
    if (geopt_inp.contains("max_force_component")) max_force_component = geopt_inp["max_force_component"];
    

    PES_optimizer::periodic_optimizer optimizer(num_atoms, init_step_size, max_history_length, minimal_step_size, subspace_tolerance);
    
    int i = 0;
    Eigen::MatrixXd forces(3, num_atoms);
    double energy;
    

    energyAndForces(mol_inp, scf_inp, energy, forces);

    Eigen::MatrixXd pos = getPositions(mol_inp);

    while (i < max_iter && forces.cwiseAbs().maxCoeff() > max_force_component) {
        std::cerr << "pos:\n";
        std::cerr << pos << + "\n";
        std::cerr << "forces:\n";
        std::cerr << forces << + "\n";
        std::cerr << i << " " << energy << " " << forces.cwiseAbs().maxCoeff() << "\n";
        optimizer.step(pos, energy, forces);
        setPositions(mol_inp, pos);
        scf_inp["initial_guess"]["type"] = "mw";
        energyAndForces(mol_inp, scf_inp, energy, forces);
        i++;
    }
    return mol_inp;

}
