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



Molecule optimize_positions(json scf_inp, json mol_inp) {

    // sqnm parameters
    Molecule mol;
    driver::init_molecule(mol_inp, mol);
    int num_atoms = mol.getNNuclei();
    double init_step_size = - 0.1;
    int max_history_length = 10;
    double minimal_step_size = 0.01;
    double subspace_tolerance = 1e-3;

    // relaxation parameters:
    int max_iter = 100;
    double max_force_component = 1e-4;

    PES_optimizer::periodic_optimizer optimizer(num_atoms, init_step_size, max_history_length, minimal_step_size, subspace_tolerance);
    
    int i = 0;
    Eigen::MatrixXd forces;
    double energy;
    

    // energyAndForces(mol_inp, scf_inp, energy, forces);

    while (i < max_iter && forces.norm() > max_force_component) {
        // optimizer.step()
        i++;
    }
    
    // return mol;

}

void energyAndForces(json mol_inp, json scf_inp, double energy, Eigen::MatrixXd &forces) {
    Molecule mol;
    driver::init_molecule(mol_inp, mol);
    auto scf_out = driver::scf::run(scf_inp, mol);
    std::cout << scf_out;
}

Eigen::MatrixXd getPositions(json mol_inp) {
    Molecule mol;
    driver::init_molecule(mol_inp, mol);
    Eigen::MatrixXd pos(3, mol.getNNuclei());

    int i = 0;
    for (const auto &coord : mol_inp["coords"]) {
        auto xyz = coord["xyz"];
        std::cout << "i would like to print";
        std::cout << xyz;
        // pos.row(i) = Eigen::Vector3d(xyz);
        i++;
    }
    return pos;
}