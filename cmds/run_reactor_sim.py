#!/usr/bin/env python3


# formal lib
import argparse
from argparse import RawTextHelpFormatter
# my lib
from microkinetic_toolkit.reactions import Reactions
#
from ase.build import fcc111
import os
import pickle


def to_kfv_dict(string):
    kfv_dict = {}
    kv_pairs = string.split(",")
    for k_tmp, v_tmp in kv_pairs:
        key = k_tmp.strip()
        va = float(v_tmp)
        kfv_dict[key] = va
    return kfv_dict


def make_parser():
    msg = """
          it runs a reactor simulation based on "
          first principle micro kinetics.
          reference: https://pubs.acs.org/doi/10.1021/acscatal.0c04104
          """
    parser = argparse.ArgumentParser(
                                description=msg,
                                fromfile_prefix_chars="@",
                                formatter_class=RawTextHelpFormatter)
    parser.add_argument("--reactions_csv", type=str, nargs="?",
                        help="csv file with reaction equation data")
    g_param_grp = parser.add_argument_group(
                                "global parameters for reactor sim")
    g_param_grp.add_argument("--T", type=float, nargs="?",
                             default=800.0, help="temperature [K]")
    g_param_grp.add_argument("--P", type=float, nargs="?",
                             default=100.0, help="total pressure [bar]")
    g_param_grp.add_argument("--p_ratio", type=to_kfv_dict, nargs="*",
                             default={"H2": 1.1, "CO2": 1.0},
                             help="partial pressure ratio\n\n"
                                  "for example, {'H2': 1.1, 'CO2': 1.0}")
    kinetic_pgrp = parser.add_argument_group(
                                "kinetic parameters")
    kinetic_pgrp.add_argument("--alpha_bep_param",
                              help="alpha in BEP relationship\n\n"
                                   "(Ea = alpha * deltaE + beta)",
                              type=float, nargs="?", default=0.6)
    kinetic_pgrp.add_argument("--beta_bep_param",
                              help="beta in BEP relationship\n\n"
                                   "(Ea = alpha * deltaE + beta)",
                              default=1.0, type=float)
    kinetic_pgrp.add_argument("--sden", type=float
                              help="site density [mol/m2]",
                              default=1.0e-4)
    kinetic_pgrp.add_argument("--v0", type=float, nargs="?", default=1.0e-3,
                              help="volumetric flow rate [m^3/sec]")
    kinetic_pgrp.add_argument("--wcat", type=float, nargs="?",
                              help="catalyst weight [mg]",
                              default=100.0e-3)
    kinetic_pgrp.add_argument("--phi", type=float, default=0.5,
                              nargs="?",
                              help="porocity")
    kinetic_pgrp.add_argument("--rho_b", type=float, nargs="?",
                              default=1.0e3, help="density")
    io_parameters = parser.add_argument_group("io paramters")
    io_parameters.add_argument("--delta_Es_file", type=str, nargs="?",
                               default="deltaEs.pickle")
    io_parameters.add_argument("--rate_graph", default=None, nargs="?",
                               type=str)
    io_parameters.add_argument("--draw_networks", default=None, nargs="?",
                               type=str)
    subparsers = parser.add_subparsers(dest="sub_opt")
    asedb_opt_parser = subparsers.add_parser("from_asedb")
    asedb_opt_parser.add_argument("--ase_db", type=str, nargs="?",
                                  default="ase.db")
    _ = subparsers.add_parser("from_emt")
    return parser


if __name__ == "__main__":
    parser = make_parser()
    args = parser.parse_args()
    # whether save files or not
    deltaEs_pickle = args.delta_Es_file
    # # read reactions from file and set Reactions
    REACTIONS_CSV = args.reactions_csv
    reactions = Reactions.from_csv(REACTIONS_CSV)
    # # define surface if you like
    reactions.surface = fcc111("Ni", size=[3, 3, 4],
                               vacuum=10.0, periodic=True)
    SUB_OPT = args.sub_opt
    if SUB_OPT == "from_asedb":
        reactions.ase_db = args.ase_db
    elif SUB_OPT == "from_emt":
        reactions.calculator = "emt"
    else:
        raise AssertionError("")
    if os.path.exists(deltaEs_pickle):
        print("loaded delta Es from pickle")
        pickle_io = open(deltaEs_pickle, "rb")
        deltaEs = pickle.load(pickle_io)
    else:
        deltaEs = reactions.get_reaction_energies()
        with open(deltaEs_pickle, "wb") as f:
            pickle.dump(deltaEs, f)
    # microkinetic analysis
    # # set parameters for kinetics
    alpha = args.alpha_bep_param
    beta = args.beta_bep_param
    # BEP relatinoship (Ea = alpha*deltaE + beta)
    bep_param = {"alpha": alpha, "beta": beta}
    sden = args.sden  # site density
    # volumetric flowrate [m^3/sec]
    v0 = args.v0
    # catalyst weight [mg]
    wcat = args.wcat
    # porocity
    phi = args.phi
    # density
    rho_b = args.rho_b
    # main function in the following.
    T = args.T
    P = args.P
    P_RATIO = args.p_ratio
    reactions.set_kinetic_parameters(
                                bep_param=bep_param,
                                sden=sden, v0=v0,
                                wcat=wcat, phi=phi, rho_b=rho_b)
    # calculate rate constant from reaction energies
    kfor = reactions.get_rate_constants(deltaEs=deltaEs, T=T)
    # # do_microkinetics
    reactions.do_microkinetics(deltaEs=deltaEs,
                               kfor=kfor, T=T, P=P,
                               ratio=P_RATIO)
    reactions.get_rate_for_graph(T=T)
    # # draw graph
    reactions.draw_network()
