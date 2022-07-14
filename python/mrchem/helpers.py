#
# MRChem, a numerical real-space code for molecular electronic structure
# calculations within the self-consistent field (SCF) approximations of quantum
# chemistry (Hartree-Fock and Density Functional Theory).
# Copyright (C) 2022 Stig Rune Jensen, Luca Frediani, Peter Wind and contributors.
#
# This file is part of MRChem.
#
# MRChem is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# MRChem is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with MRChem.  If not, see <https://www.gnu.org/licenses/>.
#
# For information on the complete list of contributors to MRChem, see:
# <https://mrchem.readthedocs.io/>
#

from .CUBEparser import write_cube_dict

# yapf: disable
SHORTHAND_FUNCTIONALS = [
    'svwn3',
    'svwn5',
    'pbe',
    'pbe0',
    'bpw91',
    'bp86',
    'b3p86',
    'b3p86-g',
    'blyp',
    'b3lyp',
    'b3lyp-g',
    'olyp',
    'kt1',
    'kt2',
    'kt3'
]
# yapf: enable
"""List of recognized shorthands for functionals"""


def write_scf_fock(user_dict, wf_dict, origin):
    fock_dict = {}

    # ZORA
    if user_dict["WaveFunction"]["relativity"].lower() == "zora":
        fock_dict["zora_operator"] = {
            "include_nuclear": user_dict["ZORA"]["include_nuclear"],
            "include_coulomb": user_dict["ZORA"]["include_coulomb"],
            "include_xc": user_dict["ZORA"]["include_xc"]
        }

    # Kinetic
    fock_dict["kinetic_operator"] = {
        "derivative": user_dict["Derivatives"]["kinetic"]
    }


    # Nuclear
    fock_dict["nuclear_operator"] = {
        "proj_prec": user_dict["Precisions"]["nuclear_prec"],
        "smooth_prec": user_dict["Precisions"]["nuclear_prec"],
        "shared_memory": user_dict["MPI"]["share_nuclear_potential"]
    }

    # Reaction
    if user_dict["WaveFunction"]["environment"].lower() == "pcm":
        fock_dict["reaction_operator"] = {
            "poisson_prec": user_dict["world_prec"],
            "kain": user_dict["PCM"]["SCRF"]["kain"],
            "max_iter": user_dict["PCM"]["SCRF"]["max_iter"],
            "optimizer": user_dict["PCM"]["SCRF"]["optimizer"],
            "dynamic_thrs": user_dict["PCM"]["SCRF"]["dynamic_thrs"],
            "density_type": user_dict["PCM"]["SCRF"]["density_type"],
            "epsilon_in": user_dict["PCM"]["Permittivity"]["epsilon_in"],
            "epsilon_out": user_dict["PCM"]["Permittivity"]["epsilon_out"],
            "formulation": user_dict["PCM"]["Permittivity"]["formulation"]
    }

    # Coulomb
    if wf_dict["method_type"] in ['hartree', 'hf', 'dft']:
        fock_dict["coulomb_operator"] = {
            "poisson_prec": user_dict["Precisions"]["poisson_prec"],
            "shared_memory": user_dict["MPI"]["share_coulomb_potential"]
        }

    # Exchange
    if wf_dict["method_type"] in ['hf', 'dft']:
        fock_dict["exchange_operator"] = {
            "poisson_prec": user_dict["Precisions"]["poisson_prec"],
            "exchange_prec": user_dict["Precisions"]["exchange_prec"]
        }

    # Exchange-Correlation
    if wf_dict["method_type"] in ['dft']:
        func_dict = []
        for line in wf_dict["dft_funcs"].split('\n'):
            sp = line.split()
            if len(sp) > 0:
                func = sp[0].lower()
                coef = [1.0]
                if len(sp) > 1:
                    coef = list(map(float, sp[1:]))
                func_dict.append({"name": func, "coef": coef[0]})
        fock_dict["xc_operator"] = {
            "shared_memory": user_dict["MPI"]["share_xc_potential"],
            "xc_functional": {
                "spin": user_dict["DFT"]["spin"],
                "cutoff": user_dict["DFT"]["density_cutoff"],
                "functionals": func_dict
            }
        }

    # External electric field
    if len(user_dict["ExternalFields"]["electric_field"]) > 0:
        fock_dict["external_operator"] = {
            "electric_field": user_dict["ExternalFields"]["electric_field"],
            "r_O": origin
        }

    return fock_dict


def write_scf_guess(user_dict, wf_dict):
    guess_str = user_dict["SCF"]["guess_type"].lower()
    guess_type = guess_str.split('_')[0]
    zeta = 0

    scf_dict = user_dict["SCF"]
    guess_prec = scf_dict["guess_prec"]
    if guess_type == 'chk':
        guess_prec = user_dict['world_prec']

    if guess_type in ['core', 'sad']:
        zeta_str = guess_str.split('_')[1]
        if zeta_str == 'sz':
            zeta = 1
        elif zeta_str == 'dz':
            zeta = 2
        elif zeta_str == 'tz':
            zeta = 3
        elif zeta_str == 'qz':
            zeta = 4
        elif zeta_str == 'gto':
            guess_type = guess_str
        else:
            print("Invalid zeta:" + zeta_str)

    file_dict = user_dict["Files"]

    if (guess_type == 'cube'):
        write_cube_dict(user_dict)

    vector_dir = file_dict["cube_vectors"]
    guess_dict = {
        "zeta": zeta,
        "prec": guess_prec,
        "type": guess_type,
        "method": wf_dict["method_name"],
        "relativity": wf_dict["relativity_name"],
        "environment": wf_dict["environment_name"],
        "external_field": wf_dict["external_name"],
        "screen": scf_dict["guess_screen"],
        "localize": scf_dict["localize"],
        "restricted": user_dict["WaveFunction"]["restricted"],
        "file_chk": scf_dict["path_checkpoint"] + "/phi_scf",
        "file_basis": file_dict["guess_basis"],
        "file_gto_p": file_dict["guess_gto_p"],
        "file_gto_a": file_dict["guess_gto_a"],
        "file_gto_b": file_dict["guess_gto_b"],
        "file_phi_p": file_dict["guess_phi_p"] + "_scf",
        "file_phi_a": file_dict["guess_phi_a"] + "_scf",
        "file_phi_b": file_dict["guess_phi_b"] + "_scf",
        "file_CUBE_p": f"{vector_dir}CUBE_p_vector.json",
        "file_CUBE_a": f"{vector_dir}CUBE_a_vector.json",
        "file_CUBE_b": f"{vector_dir}CUBE_b_vector.json"

    }
    return guess_dict


def write_scf_solver(user_dict, wf_dict):
    # SCF precisions and thresholds
    start_prec = user_dict["SCF"]["start_prec"]
    final_prec = user_dict["SCF"]["final_prec"]
    if final_prec < 0.0:
        final_prec = user_dict["world_prec"]
    if start_prec < 0.0:
        start_prec = final_prec

    scf_dict = user_dict["SCF"]
    solver_dict = {
        "method": wf_dict["method_name"],
        "relativity": wf_dict["relativity_name"],
        "environment": wf_dict["environment_name"],
        "external_field": wf_dict["external_name"],
        "kain": scf_dict["kain"],
        "max_iter": scf_dict["max_iter"],
        "rotation": scf_dict["rotation"],
        "localize": scf_dict["localize"],
        "file_chk": scf_dict["path_checkpoint"] + "/phi_scf",
        "checkpoint": scf_dict["write_checkpoint"],
        "start_prec": start_prec,
        "final_prec": final_prec,
        "energy_thrs": scf_dict["energy_thrs"],
        "orbital_thrs": scf_dict["orbital_thrs"],
        "helmholtz_prec": user_dict["Precisions"]["helmholtz_prec"]
    }

    return solver_dict


def write_scf_properties(user_dict, origin):
    prop_dict = {}
    if user_dict["Properties"]["dipole_moment"]:
        prop_dict["dipole_moment"] = {}
        prop_dict["dipole_moment"]['dip-1'] = {
            "operator": "h_e_dip",
            "precision": user_dict["world_prec"],
            "r_O": origin
        }
    if user_dict["Properties"]["quadrupole_moment"]:
        prop_dict["quadrupole_moment"] = {}
        prop_dict["quadrupole_moment"]['quad-1'] = {
            "operator": "h_e_quad",
            "precision": user_dict["world_prec"],
            "r_O": origin
        }
    if user_dict["Properties"]["geometric_derivative"]:
        prop_dict["geometric_derivative"] = {}
        prop_dict["geometric_derivative"]["geom-1"] = {
            "operator": "h_nuc_grad",
            "precision": user_dict["world_prec"],
            "smoothing": user_dict["Precisions"]["nuclear_prec"],
        }
    return prop_dict


def write_scf_plot(user_dict):
    plot_dict = {}
    if user_dict["Properties"]["plot_density"] or len(user_dict["Properties"]["plot_orbitals"]):
        plot_dict["orbitals"] = user_dict["Properties"]["plot_orbitals"]
        plot_dict["density"] = user_dict["Properties"]["plot_density"]
        plot_dict["plotter"] = user_dict["Plotter"]
        if user_dict["world_unit"] == "angstrom":
            plot_dict["plotter"] = {
                k: [user_dict['Constants']['angstrom2bohrs'] * r for r in plot_dict["plotter"][k]]
                for k in plot_dict["plotter"].keys()
            }
    return plot_dict


def write_rsp_calc(omega, user_dict, origin):
    wf_dict = parse_wf_method(user_dict)
    if not wf_dict["relativity_name"] in ["None", "Off"]:
        raise RuntimeError("Linear response not available: " + wf_dict["relativity_name"])

    rsp_dict = user_dict["Response"]
    file_dict = user_dict["Files"]

    rsp_calc = {}
    rsp_calc["frequency"] = omega
    rsp_calc["dynamic"] = (omega > 1.0e-12)
    rsp_calc["fock_operator"] = write_rsp_fock(user_dict, wf_dict)
    rsp_calc["unperturbed"] = {
        "precision": user_dict["world_prec"],
        "localize": rsp_dict["localize"],
        "fock_operator": write_scf_fock(user_dict, wf_dict, origin)
    }

    guess_str = rsp_dict["guess_type"].lower()
    guess_type = guess_str.split('_')[0]
    guess_prec = rsp_dict["guess_prec"]
    if guess_type == 'chk':
        guess_prec = user_dict['world_prec']

    rsp_calc["components"] = []
    for d in [0, 1, 2]:
        rsp_comp = {}
        rsp_comp["initial_guess"] = {
            "prec": guess_prec,
            "type": guess_type,
            "file_chk_x": rsp_dict["path_checkpoint"] + "/X_rsp_" + str(d),
            "file_chk_y": rsp_dict["path_checkpoint"] + "/Y_rsp_" + str(d),
            "file_x_p": file_dict["guess_x_p"] + "_rsp_" + str(d),
            "file_x_a": file_dict["guess_x_a"] + "_rsp_" + str(d),
            "file_x_b": file_dict["guess_x_b"] + "_rsp_" + str(d),
            "file_y_p": file_dict["guess_y_p"] + "_rsp_" + str(d),
            "file_y_a": file_dict["guess_y_a"] + "_rsp_" + str(d),
            "file_y_b": file_dict["guess_y_b"] + "_rsp_" + str(d)
        }
        if rsp_dict["write_orbitals"]:
            path_orbitals = rsp_dict["path_orbitals"]
            rsp_comp["write_orbitals"] = {
                "file_x_p": path_orbitals + "/X_p_rsp_" + str(d),
                "file_x_a": path_orbitals + "/X_a_rsp_" + str(d),
                "file_x_b": path_orbitals + "/X_b_rsp_" + str(d),
                "file_y_p": path_orbitals + "/Y_p_rsp_" + str(d),
                "file_y_a": path_orbitals + "/Y_a_rsp_" + str(d),
                "file_y_b": path_orbitals + "/Y_b_rsp_" + str(d)
            }
        if rsp_dict["run"][d]:
            rsp_comp["rsp_solver"] = write_rsp_solver(user_dict, wf_dict, d)
        rsp_calc["components"].append(rsp_comp)
    return rsp_calc


def write_rsp_fock(user_dict, wf_dict):
    fock_dict = {}

    # Coulomb
    if wf_dict["method_type"] in ['hartree', 'hf', 'dft']:
        fock_dict["coulomb_operator"] = {
            "poisson_prec": user_dict["Precisions"]["poisson_prec"],
            "shared_memory": user_dict["MPI"]["share_coulomb_potential"]
        }

    # Exchange
    if wf_dict["method_type"] in ['hf', 'dft']:
        fock_dict["exchange_operator"] = {
            "poisson_prec": user_dict["Precisions"]["poisson_prec"],
            "exchange_prec": user_dict["Precisions"]["exchange_prec"]
        }

    # Exchange-Correlation
    if wf_dict["method_type"] in ['dft']:
        func_dict = []
        for line in wf_dict["dft_funcs"].split('\n'):
            sp = line.split()
            if len(sp) > 0:
                func = sp[0].lower()
                coef = [1.0]
                if len(sp) > 1:
                    coef = list(map(float, sp[1:]))
                func_dict.append({"name": func, "coef": coef[0]})
        fock_dict["xc_operator"] = {
            "shared_memory": user_dict["MPI"]["share_xc_potential"],
            "xc_functional": {
                "spin": user_dict["DFT"]["spin"],
                "cutoff": user_dict["DFT"]["density_cutoff"],
                "functionals": func_dict
            }
        }

    return fock_dict


def write_rsp_solver(user_dict, wf_dict, d):
    # Response precisions and thresholds
    start_prec = user_dict["Response"]["start_prec"]
    final_prec = user_dict["Response"]["final_prec"]
    if final_prec < 0.0:
        final_prec = user_dict["world_prec"]
    if start_prec < 0.0:
        start_prec = final_prec

    rsp_dict = user_dict["Response"]
    solver_dict = {
        "method": wf_dict["method_name"],
        "kain": rsp_dict["kain"],
        "max_iter": rsp_dict["max_iter"],
        "file_chk_x": rsp_dict["path_checkpoint"] + "/X_rsp_" + str(d),
        "file_chk_y": rsp_dict["path_checkpoint"] + "/Y_rsp_" + str(d),
        "checkpoint": rsp_dict["write_checkpoint"],
        "start_prec": start_prec,
        "final_prec": final_prec,
        "orbital_thrs": user_dict["Response"]["orbital_thrs"],
        "property_thrs": user_dict["Response"]["property_thrs"],
        "helmholtz_prec": user_dict["Precisions"]["helmholtz_prec"],
        "orth_prec": 1.0e-14
    }
    return solver_dict


def parse_wf_method(user_dict):
    method_name = ''
    restricted = user_dict["WaveFunction"]["restricted"]
    method_type = user_dict["WaveFunction"]["method"].lower()
    dft_funcs = user_dict["DFT"]["functionals"].lower()
    if method_type in ['core']:
        method_name = 'Core Hamiltonian'
    elif method_type in ['hartree']:
        method_name = 'Hartree'
    elif method_type in ['hf', 'hartree-fock', 'hartreefock']:
        method_name = 'Hartree-Fock'
        method_type = 'hf'
    elif method_type in ['dft']:
        method_name = 'DFT'
    elif method_type in ['lda']:
        method_name = 'DFT (SVWN5)'
        dft_funcs = 'svwn5'
        method_type = 'dft'
    elif method_type in SHORTHAND_FUNCTIONALS:
        method_name = 'DFT (' + method_type.upper() + ')'
        dft_funcs = method_type
        method_type = 'dft'
    else:
        raise RuntimeError(
            f"Invalid wavefunction method {user_dict['WaveFunction']['method']}"
        )

    # Determine relativity name label for print outs to the output file
    relativity_name = "None"
    if user_dict["WaveFunction"]["relativity"].lower() in ["none"]:
        user_dict["WaveFunction"]["relativity"] = "off"
        user_dict["ZORA"]["include_nuclear"] = False
        user_dict["ZORA"]["include_coulomb"] = False
        user_dict["ZORA"]["include_xc"] = False

    if user_dict["WaveFunction"]["relativity"].lower() in ["nzora"]:
        user_dict["WaveFunction"]["relativity"] = "zora"
        user_dict["ZORA"]["include_nuclear"] = True
        user_dict["ZORA"]["include_coulomb"] = False
        user_dict["ZORA"]["include_xc"] = False

    if user_dict["WaveFunction"]["relativity"].lower() in ["zora"]:
        components = [
            user_dict["ZORA"]["include_nuclear"],
            user_dict["ZORA"]["include_coulomb"],
            user_dict["ZORA"]["include_xc"]
        ]
        names = [
            "V_nuc",
            "J",
            "V_xc"
        ]

        if any(components):
            zora_terms = " + ".join([name for name, comp in zip(names, components) if comp])
            relativity_name = "ZORA (" + zora_terms + ")"
        else:
            raise RuntimeError("ZORA selected, but no ZORA potential included")

        if user_dict["ZORA"]["include_xc"] and not restricted:
            raise RuntimeError("ZORA (V_xc) not available for unrestricted wavefunctions")


    # Determine environment name label for print outs to the output file
    environment_name = "None"
    if user_dict["WaveFunction"]["environment"].lower() == "pcm":
        environment_name = "PCM"


    # Determine external name label for print outs to the output file
    ext_dict = user_dict['ExternalFields']
    has_external_fields = len(ext_dict['electric_field']) > 0

    external_name = "None"
    if has_external_fields:
        # If no external fields, then the list will be empty
        # Need to catch the exception and store placeholders
        try:
            x, y, z = ext_dict['electric_field']
        except ValueError:
            x, y, z = None, None, None  # Useless placeholders

        # Labels to aggregate
        external_name = f'Electric field ({x}, {y}, {z})'


    wf_dict = {
        "relativity_name": relativity_name,
        "environment_name": environment_name,
        "external_name": external_name,
        "method_name": method_name,
        "method_type": method_type,
        "dft_funcs": dft_funcs
    }
    return wf_dict
