#!/usr/bin/env python
from __future__ import division

# Interface pipeline: designs binders with HIS that h-bond the target, scored by pH sensitivity.
# Needs the hacked ProteinMPNN (set RF_DIFFUSION_PATH).
#
# Usage: ./his_ph_interface_design.py pdb1.pdb pdb2.pdb [options]
#    or: ./his_ph_interface_design.py -in:file:silent my.silent [options]


import os
import sys
import math

import os
import sys

import pyrosetta
from pyrosetta import *
from pyrosetta.rosetta import *

import functools

# ProteinMPNN is run in one of two ways:
#  1. Stock ProteinMPNN as an external program. It is found through PROTEIN_MPNN_PATH (a clone of
#      https://github.com/dauparas/ProteinMPNN) or a ProteinMPNN folder next to this script.
#      Optionally set PROTEIN_MPNN_PYTHON to a python that has torch. See his_ph_interface_mpnn.py
#  2. The author's hacked ProteinMPNN (inference/bcov_hacks/hacked_protein_mpnn_run.py in an rf_diffusion checkout), in-process.
#      Set RF_DIFFUSION_PATH. It needs torch and openfold.
# Stock ProteinMPNN wins if both are available.
import his_ph_interface_mpnn
hacked_protein_mpnn_run = None
if his_ph_interface_mpnn.find_mpnn_path() is None and "-h" not in sys.argv and "--help" not in sys.argv:
    if not os.environ.get("RF_DIFFUSION_PATH"):
        sys.exit("Can't find ProteinMPNN. Set PROTEIN_MPNN_PATH to a clone of https://github.com/dauparas/ProteinMPNN")
    sys.path.append(os.environ["RF_DIFFUSION_PATH"])
    import inference.bcov_hacks.hacked_protein_mpnn_run as hacked_protein_mpnn_run

import importlib





import npose_util_pyrosetta as nup
import npose_util as nu

import numpy as np
from collections import defaultdict
import time
import argparse
import itertools
import subprocess
import time
import re

import scipy.stats
import shlex

import pandas as pd

def _hbedge_from_lowmem(hb_graph, lowmem_edge, node_ind):
    # The edge list iterates over LowMemEdges. We need the full HBondEdge (with the hbonds).
    if hasattr(hb_graph, 'HBondEdge_from_LowMemEdge'):  # pyrosetta with the hbond graph patch
        return hb_graph.HBondEdge_from_LowMemEdge(lowmem_edge)
    # Published pyrosetta: look the full edge up from its two node indices
    return hb_graph.find_edge(node_ind, lowmem_edge.get_other_ind(node_ind))

# import pyRMSD.RMSDCalculator

init("-beta_nov16 -in:file:silent_struct_type binary -keep_input_scores false -mute all"
    )


# setup pH mode
basic.options.set_boolean_option('pH:pH_mode', True)
pose = pose_from_sequence('H')
scorefxn = get_fa_scorefxn()
protocols.toolbox.pose_manipulation.repack_this_residue(1, pose, scorefxn)
basic.options.set_boolean_option('pH:pH_mode', False)



parser = argparse.ArgumentParser()


parser.add_argument('-his_in_groups_of', type=int, default=3)
parser.add_argument('-first_round_n_sets', type=int, default=10)
parser.add_argument('-second_round_n_sets', type=int, default=10)
parser.add_argument('-min_his_score', type=int, default=5)


parser.add_argument("-in:file:silent", type=str, default="")
parser.add_argument("pdbs", type=str, nargs="*")
# parser.add_argument("-num_to_output", type=int, default=5)
parser.add_argument("-num_per_input", type=int, default=20)
parser.add_argument("-mpnn_seqs", type=int, default=1)
parser.add_argument("-mpnn_temps", type=str, default="0.001 0.01 0.1")
# parser.add_argument("-mpnn_n_cycles", type=int, default=1)
# parser.add_argument("-mpnn_single_chain", action='store_true')


parser.add_argument("-longest_loop", type=int, default=5)
parser.add_argument("-percent_core_scn", type=float, default=0.22)
parser.add_argument("-micro_helices", type=int, default=0)
parser.add_argument("-any_core_9", type=float, default=1)
parser.add_argument("-other_hits_9", type=float, default=1)

parser.add_argument("-do_monomer_filter", action="store_true")
parser.add_argument("-crappy_interface_mode", action="store_true")
parser.add_argument("-trim_distance", type=float, default=15)

parser.add_argument("-sap_limit", type=float, default=100)
parser.add_argument("-filter_on_sap_too", action="store_true")

parser.add_argument("-dont_filter_patch", action="store_true")

parser.add_argument("-neg_netc_bias_level", type=float, default=0.0)
parser.add_argument("-surface_ala_weight", type=float, default=-0.75)
parser.add_argument("-interface_desap_bcov", type=float, default=0)
parser.add_argument("-interface_desap_sap", type=float, default=0)

parser.add_argument("-force_seq_by_tag", type=str, default="")

parser.add_argument("-pssm_by_tag", type=str, default="") # just a dictionary with seqpos then letter then value
parser.add_argument("-pssm_multiplier", type=float, default=1)

parser.add_argument("-hbnet_lock_identities", action='store_true')
parser.add_argument("-hbnet_lock_identities_strict", action='store_true')

parser.add_argument("-swap_target", type=str, default="")

parser.add_argument("-patchdock_res_regions_weights", type=str, 
    default="", 
    help="1,2,3:0.5;4,5,6:2.3;100,101,102:3.14")


parser.add_argument("-tied_target", type=str, default='')


# assert(False fix surface ala, dont pack loops, also more memory, clean exit if none pas)

args = parser.parse_args(sys.argv[1:])

pdbs = args.pdbs
silent = args.__getattribute__("in:file:silent")




scorefxn = get_fa_scorefxn()
scorefxn_fa_atr = core.scoring.ScoreFunctionFactory.create_score_function("none")
scorefxn_fa_atr.set_weight(core.scoring.fa_atr, 1)

scorefxn_none = core.scoring.ScoreFunctionFactory.create_score_function("none")


chainA = core.select.residue_selector.ChainSelector("A")
chainB = core.select.residue_selector.ChainSelector("B")
interface_on_A = core.select.residue_selector.NeighborhoodResidueSelector(chainB, 10.0, False)
interface_on_B = core.select.residue_selector.NeighborhoodResidueSelector(chainA, 10.0, False)
big_interface_on_A = core.select.residue_selector.NeighborhoodResidueSelector(chainB, 14.0, False)
big_interface_on_B = core.select.residue_selector.NeighborhoodResidueSelector(chainA, 14.0, False)
interface_by_vector = core.select.residue_selector.InterGroupInterfaceByVectorSelector(interface_on_A, interface_on_B)
interface_by_vector.cb_dist_cut(11)
interface_by_vector.cb_dist_cut(5.5)
interface_by_vector.vector_angle_cut(75)
interface_by_vector.vector_dist_cut(9)

new_target = None
if ( args.swap_target ):
    new_target = pose_from_file(args.swap_target)

tied_target = None
if args.tied_target:
    tied_target = pose_from_file(args.tied_target)


A_or_big_B = core.select.residue_selector.OrResidueSelector(chainA, big_interface_on_B)


force_seq_by_tag = {}
if ( args.force_seq_by_tag != "" ):
    with open(args.force_seq_by_tag) as f:
        for line in f:
            line = line.strip()
            if (len(line) == 0 ):
                continue
            sp = line.split()
            assert(len(sp) == 2)
            force_seq_by_tag[sp[1]] = sp[0]

    assert(len(force_seq_by_tag) > 0)

def get_mpnn(design_chains='A'):
    if hacked_protein_mpnn_run is None:
        return functools.partial(his_ph_interface_mpnn.run_stock_mpnn, design_chains)

    argparser = hacked_protein_mpnn_run.load_argparser()
    args = argparser.parse_args(shlex.split(f'--pdb_path_chains="{design_chains}" --out_folder ./ --path_to_model_weights= --omit_AAs C'))

    mpnn = functools.partial(hacked_protein_mpnn_run.main, args)
    return mpnn


design_chains = 'A' if tied_target is None else 'A C'
mpnn = get_mpnn(design_chains=design_chains)


surface_bias_AA_dict = {}
if ( args.neg_netc_bias_level != 0 ):
    surface_bias_AA_dict['E'] = args.neg_netc_bias_level
    surface_bias_AA_dict['R'] = -args.neg_netc_bias_level/2
    surface_bias_AA_dict['K'] = -args.neg_netc_bias_level/2
surface_bias_AA_dict['A'] = args.surface_ala_weight



sap_hydrophobicity_values = {
    "A": 0.116,
    "C": 0.18,
    "D": -0.472,
    "E": -0.457,
    "F": 0.5,
    "G": 0.001,
    "H": -0.335,
    "I": 0.443,
    "K": -0.217,
    "L": 0.443,
    "M": 0.238,
    "N": -0.264,
    "P": 0.211,
    "Q": -0.249,
    "R": 0.0, # R modified to 0 from -0.5 for design
    "S": -0.141,
    "T": -0.05,
    "V": 0.325,
    "W": 0.378,
    "Y": 0.38,
}

arbitrary_bcov_desap_weights = {
'D':0.5,
'E':1,
'H':1,
'K':1,
'N':0.5,
'Q':1,
'R':1,
'S':0.25,
'T':0.25,
'Y':0.1,
}



interface_bias_AA_dict = {}

if args.interface_desap_sap != 0:
    for letter in sap_hydrophobicity_values:
        to_add = args.interface_desap_sap * sap_hydrophobicity_values[letter] * -1 # the values in the table are backwards
        if letter in interface_bias_AA_dict:
            interface_bias_AA_dict[letter] += to_add
        else:
            interface_bias_AA_dict[letter] = to_add

if args.interface_desap_bcov != 0:
    for letter in arbitrary_bcov_desap_weights:
        to_add = args.interface_desap_bcov * arbitrary_bcov_desap_weights[letter]
        if letter in interface_bias_AA_dict:
            interface_bias_AA_dict[letter] += to_add
        else:
            interface_bias_AA_dict[letter] = to_add







pssms = None
if args.pssm_by_tag != "":
    pssms = {}
    with open(args.pssm_by_tag) as f:
        for line in f:
            if len(line) == 0:
                continue
            sp = line.split()
            assert(len(sp) == 2)
            pssms[sp[1]] = sp[0].replace("_", " ")


all_patch_indices = set()

selectors = []
filters = []
names = []
cp_weights = []

if args.patchdock_res_regions_weights == '':
    args.dont_filter_patch = True
    args.patchdock_res_regions_weights = '1:1'


selectors.append('<Chain name="chainA" chains="A"/>')
selectors.append('<Chain name="chainB" chains="B"/>')

for i, res_and_weight in enumerate(args.patchdock_res_regions_weights.split(";")):
    if ( len(res_and_weight) == 0 ):
        continue
    res_string, weight_string = res_and_weight.split(":")

    cp_weights.append( float(weight_string))

    name = "contact_patch%i"%i

    selectors.append(f'<Slice name="{name}" indices="{res_string}" selector="chainB" />')
    filters.append(f'<ContactMolecularSurface name="{name}" distance_weight="0.5" target_selector="{name}"'
                                                                    +' apolar_target="true" binder_selector="chainA" confidence="0" />')
    names.append(name)

    for idx in res_string.split(","):
        all_patch_indices.add(int(idx))

selectors_string = "\n".join(selectors)
filters_string = "\n".join(filters)
xml = f'''
<RESIDUE_SELECTORS>
{selectors_string}
</RESIDUE_SELECTORS>
<FILTERS>
{filters_string}
</FILTERS>
'''

objs = protocols.rosetta_scripts.XmlObjects.create_from_string(xml)

cp_filters = []
for name in names:
    cp_filters.append(objs.get_filter(name))

def calc_cp_score(pose):
    score = 0
    for cp_filter, weight in zip(cp_filters, cp_weights):
        score += weight * cp_filter.report_sm(pose)

    return score


    
def fix_scorefxn(sfxn, allow_double_bb=False):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    sfxn.set_energy_method_options(opts)

# fix_scorefxn(scorefxn)


def my_rstrip(string, strip):
    if (string.endswith(strip)):
        return string[:-len(strip)]
    return string


def add_to_score_map(og_map, to_add, prefix, suffix=""):
    for name, score in list(to_add.items()):    # this iterator is broken. use list()
        og_map[prefix + name + suffix] = score

def move_chainA_far_away(pose):
    pose = pose.clone()
    sel = core.select.residue_selector.ChainSelector("A")
    subset = sel.apply(pose)

    x_unit = numeric.xyzVector_double_t(1, 0, 0)
    far_away = numeric.xyzVector_double_t(10000, 0, 0)

    protocols.toolbox.pose_manipulation.rigid_body_move(x_unit, 0, far_away, pose, subset)

    return pose


def which_chain(pose, resi):
    for i in range(1, pose.num_chains()+1):
        if ( pose.conformation().chain_begin(i) <= resi and pose.conformation().chain_end(i) >= resi):
            return i
    assert(False)

def get_monomer_score(pose, scorefxn):
    pose = pose.split_by_chain()[1]
    return scorefxn(pose)


def get_filter_by_name(filtername):
    the_filter = objs.get_filter(filtername)

    # Get rid of stochastic filter
    if ( isinstance(the_filter, pyrosetta.rosetta.protocols.filters.StochasticFilter) ):
        the_filter = the_filter.subfilter()

    return the_filter

def add_filter_to_results(pose, filtername, out_score_map):
    filter = get_filter_by_name(filtername)
    print("protocols.rosetta_scripts.ParsedProtocol.REPORT: ============Begin report for " + filtername + "=======================")
    if (isinstance(filter, protocols.simple_filters.ShapeComplementarityFilter)):
        value = filter.compute(pose)
        out_score_map[filtername] = value.sc
        out_score_map[filtername+"_median_dist"] = value.distance
    else:
        value = filter.report_sm(pose)
        out_score_map[filtername] = value
    print("============End report for " + filtername + "=======================")

def score_with_these_filters(pose, filters, out_score_map):
    for filtername in filters:
        add_filter_to_results(pose, filtername, out_score_map)

def cmd(command, wait=True):
    # print ""
    # print command
    the_command = subprocess.Popen(command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True)
    if (not wait):
        return
    the_stuff = the_command.communicate()
    return str(the_stuff[0]) + str(the_stuff[1])

def atid(resnum, atno):
    return core.id.AtomID( atno, resnum )

abego_man = core.sequence.ABEGOManager()
def get_abego(pose, seqpos):
    return abego_man.index2symbol(abego_man.torsion2index_level1( pose.phi(seqpos), pose.psi(seqpos), pose.omega(seqpos)))

def dump_region(pose, start, end, name):
    residues = utility.vector1_unsigned_long()
    for i in range(start, end ):
        residues.append(i)

    to_dump = core.pose.Pose()
    core.pose.pdbslice(to_dump, pose, residues)
    pdbinfo = core.pose.PDBInfo( to_dump )
    to_dump.pdb_info(pdbinfo)
    to_dump.dump_pdb(name)


def get_per_atom_sasa(pose):
    atoms = core.id.AtomID_Map_bool_t()
    atoms.resize(pose.size())
    for i in range(1, pose.size()+1):
        atoms.resize( i, pose.residue(i).natoms(), True)
    surf_vol = core.scoring.packing.get_surf_vol( pose, atoms, 2.8)
    # print(surf_vol.tot_surf)
    # print(surf_vol.surf(2, 1))  # this is per atom sasa (residue 2, atom 1)
    return surf_vol

# this is 1 indexed with the start and end as xx
# and HHHHHH turns identified
def better_dssp(pose, length=-1):
    if ( length < 0 ):
        length = pose.size()

    dssp = core.scoring.dssp.Dssp(pose)
    dssp.dssp_reduced()
    the_dssp = "x" + dssp.get_dssp_secstruct()
    the_dssp = list(the_dssp)
    the_dssp[1] = "x"
    the_dssp[-1] = "x"
    the_dssp[2] = "x"
    the_dssp[-2] = "x"
    the_dssp[3] = "x"
    the_dssp[-3] = "x"
    the_dssp = "".join(the_dssp)

    my_dssp = "x"

    for seqpos in range(1, length+1):
        abego = get_abego(pose, seqpos)
        this_dssp = the_dssp[seqpos]
        if ( the_dssp[seqpos] == "H" and abego != "A" ):
            # print("!!!!!!!!!! Dssp - abego mismatch: %i %s %s !!!!!!!!!!!!!!!"%(seqpos, the_dssp[seqpos], abego))

            # This is the Helix-turn-helix HHHH case. See the test_scaffs folder
            if ( abego == "B" ):
                this_dssp = "L"

        my_dssp += this_dssp

    return my_dssp


def get_consensus(letters):
    counts = defaultdict(lambda : 0, {})
    for letter in letters:
        counts[letter] += 1

    maxx_letter = 0
    maxx = 0
    for key in counts:
        if ( counts[key] > maxx ):
            maxx = counts[key]
            maxx_letter = key
    return maxx_letter

# this is 1 indexed with the start and end with loops converted to nearby dssp
# and HHHHHH turns identified
def better_dssp2(pose, length=-1):
    if ( length < 0 ):
        length = pose.size()

    dssp = core.scoring.dssp.Dssp(pose)
    dssp.dssp_reduced()
    the_dssp = "x" + dssp.get_dssp_secstruct()
    the_dssp = list(the_dssp)

    n_consensus = get_consensus(the_dssp[3:6])

    the_dssp[1] = n_consensus
    the_dssp[2] = n_consensus
    the_dssp[3] = n_consensus
    the_dssp[4] = n_consensus
    the_dssp[5] = n_consensus

    c_consensus = get_consensus(the_dssp[-5:-2])

    the_dssp[-1] = c_consensus
    the_dssp[-2] = c_consensus
    the_dssp[-3] = c_consensus
    the_dssp[-4] = c_consensus
    the_dssp[-5] = c_consensus

    the_dssp = "".join(the_dssp)

    my_dssp = "x"

    for seqpos in range(1, length+1):
        abego = get_abego(pose, seqpos)
        this_dssp = the_dssp[seqpos]
        if ( the_dssp[seqpos] == "H" and abego != "A" ):
            # print("!!!!!!!!!! Dssp - abego mismatch: %i %s %s !!!!!!!!!!!!!!!"%(seqpos, the_dssp[seqpos], abego))

            # This is the Helix-turn-helix HHHH case. See the test_scaffs folder
            if ( abego == "B" or abego == "E" and seqpos > 5 and seqpos < len(the_dssp)-5 ):
                this_dssp = "L"

        my_dssp += this_dssp

    return my_dssp


# this is 1 indexed with the start and end with loops converted to nearby dssp
# and HHHHHH turns identified
def better_dssp3(pose, length=-1, force_consensus=None, consensus_size=6):
    if ( length < 0 ):
        length = pose.size()

    dssp = core.scoring.dssp.Dssp(pose)
    dssp.dssp_reduced()
    the_dssp = "x" + dssp.get_dssp_secstruct()[:length]
    the_dssp = list(the_dssp)

    n_consensus = get_consensus(the_dssp[3:consensus_size+1])
    if ( not force_consensus is None ):
        n_consensus = force_consensus

    for i in range(1, consensus_size+1):
        the_dssp[i] = n_consensus

    c_consensus = get_consensus(the_dssp[-(consensus_size):-2])
    if ( not force_consensus is None ):
        c_consensus = force_consensus

    for i in range(1, consensus_size+1):
        the_dssp[-i] = c_consensus

    the_dssp = "".join(the_dssp)

    # print(the_dssp)

    my_dssp = "x"

    for seqpos in range(1, length+1):
        abego = get_abego(pose, seqpos)
        this_dssp = the_dssp[seqpos]
        if ( the_dssp[seqpos] == "H" and abego != "A" ):
            # print("!!!!!!!!!! Dssp - abego mismatch: %i %s %s !!!!!!!!!!!!!!!"%(seqpos, the_dssp[seqpos], abego))

            # This is the Helix-turn-helix HHHH case. See the test_scaffs folder
            if ( (abego == "B" or abego == "E") and seqpos > consensus_size and seqpos < len(the_dssp)-consensus_size ):
                this_dssp = "L"

        my_dssp += this_dssp

    # print(my_dssp)

    return my_dssp



# assumes dssp starts with X and removes it
def get_ss_elements2(dssp):
    assert(dssp[0] == "x")
    ss_elements = []

    offset = 0
    ilabel = -1
    for label, group in itertools.groupby(dssp):
        ilabel += 1
        this_len = sum(1 for _ in group)
        next_offset = offset + this_len

        ss_elements.append( (label, offset, next_offset-1))

        offset = next_offset
    return ss_elements[1:]



def superposition_transform_with_weight(input_coords_move, input_coords_to, rotation_matrix, move_com, ref_com, atom_weight):

    coords_move = utility.vector1_numeric_xyzVector_double_t()
    coords_to = utility.vector1_numeric_xyzVector_double_t()

    assert(len(input_coords_move) == len(input_coords_to))
    assert(len(input_coords_move) == len(atom_weight))

    for i in range(len(input_coords_move)):
        for j in range(atom_weight[i]):
            coords_move.append(input_coords_move[i+1])
            coords_to.append(input_coords_to[i+1])

    protocols.toolbox.superposition_transform( coords_move, coords_to, rotation_matrix, move_com, ref_com )


# align two things using sequence and accepting gaps
def pymol_align( move_pose, to_pose, sel_move=None, sel_to=None, atoms=["N", "CA", "C"], throw_away=0.1, extend_penalty=-1,
        to_pose_upweight_mask=None, move_pose_upweight_mask=None, return_all_move_to=False ):

    if ( to_pose_upweight_mask is None ):
        to_pose_upweight_mask = np.zeros(to_pose.size()+1, bool)

    if ( move_pose_upweight_mask is None ):
        move_pose_upweight_mask = np.zeros(move_pose.size()+1, bool)

    if ( not sel_move is None ):
        move_res = np.array(list(core.select.get_residues_from_subset(sel_move)))
    else:
        move_res = np.array(list(range(1, move_pose.size()+1)))

    if ( not sel_to is None ):
        to_res = np.array(list(core.select.get_residues_from_subset(sel_to)))
    else:
        to_res = np.array(list(range(1, to_pose.size()+1)))

    seq_move = "x" + move_pose.sequence()
    seq_to = "x" + to_pose.sequence()

    seq_move = "".join(np.array(list(seq_move))[move_res])
    seq_to = "".join(np.array(list(seq_to))[to_res])

    from Bio import pairwise2
    alignment = align_move, align_to, idk1, idk2, idk3 = pairwise2.align.globalxs(seq_move,seq_to, -2, extend_penalty, penalize_end_gaps=False)[0]
    # print(align_move, align_to)
    # print(alignment)

    all_move_to = []
    move_to_pairs = []
    coords_move = utility.vector1_numeric_xyzVector_double_t()
    coords_to = utility.vector1_numeric_xyzVector_double_t()
    atom_weight = []

    i_move = 0
    i_to = 0
    for i in range(len(align_move)):
        if ( align_move[i] == align_to[i] ):

            seqpos_move = move_res[i_move]
            seqpos_to = to_res[i_to]

            move_to_pairs.append((seqpos_move, seqpos_to))
            all_move_to.append((seqpos_move, seqpos_to))


            weight = 1
            if ( to_pose_upweight_mask[seqpos_to] ):
                weight = 5
            if ( move_pose_upweight_mask[seqpos_move] ):
                weight = 5


            for atom in atoms:
                coords_move.append(move_pose.residue(seqpos_move).xyz(atom))
                coords_to.append(to_pose.residue(seqpos_to).xyz(atom))
                atom_weight.append(weight)
        else:

            if ( align_move[i] != "-" and align_to[i] != "-" ):

                seqpos_move = move_res[i_move]
                seqpos_to = to_res[i_to]

                all_move_to.append((seqpos_move, seqpos_to))



        if ( align_move[i] != "-" ):
            i_move += 1
        if ( align_to[i] != "-" ):
            i_to += 1

    move_pose_copy = move_pose.clone()


    rmsd = 0

    distances = []

    if ( len(move_to_pairs) > 0 ):

        rotation_matrix = numeric.xyzMatrix_double_t()
        move_com = numeric.xyzVector_double_t()
        ref_com = numeric.xyzVector_double_t()

        superposition_transform_with_weight(coords_move, coords_to, rotation_matrix, move_com, ref_com, atom_weight)
        # protocols.toolbox.superposition_transform( coords_move, coords_to, rotation_matrix, move_com, ref_com )

        protocols.toolbox.apply_superposition_transform(move_pose, rotation_matrix, move_com, ref_com)

        for seqpos_move, seqpos_to in move_to_pairs:
            for atom in atoms:
                distance = move_pose.residue(seqpos_move).xyz(atom).distance_squared(to_pose.residue(seqpos_to).xyz(atom))
                rmsd += distance
                distances.append(distance)

        rmsd /= len(move_to_pairs)*len(atoms)
        rmsd = np.sqrt(rmsd)

    move_pose = move_pose_copy

    distances = np.array(distances)

    # print("Initial RMSD: %.3f"%rmsd)

    cutoff = np.percentile(distances, 100 - throw_away * 10)
    # print("Cutoff %.3f"%cutoff)

    mask = distances <= cutoff

    # print(mask.sum(), len(mask))

    coords_move_old = list(coords_move)
    coords_to_old = list(coords_to)
    atom_weight_old = list(atom_weight)
    # move_to_pairs_old = move_to_pairs

    # move_to_pairs = []
    coords_move = utility.vector1_numeric_xyzVector_double_t()
    coords_to = utility.vector1_numeric_xyzVector_double_t()
    atom_weight = []

    for i in range(len(coords_move_old)):
        if ( not mask[i] ):
            continue
        coords_move.append(coords_move_old[i])
        coords_to.append(coords_to_old[i])
        atom_weight.append(atom_weight_old[i])
        # move_to_pairs.append(move_to_pairs_old[i])

    # print(len(coords_move), len(coords_move_old))

    rmsd = 0
    imask = -1
    if ( len(move_to_pairs) > 0 ):

        rotation_matrix = numeric.xyzMatrix_double_t()
        move_com = numeric.xyzVector_double_t()
        ref_com = numeric.xyzVector_double_t()

        superposition_transform_with_weight( coords_move, coords_to, rotation_matrix, move_com, ref_com, atom_weight )

        protocols.toolbox.apply_superposition_transform(move_pose, rotation_matrix, move_com, ref_com)

        for seqpos_move, seqpos_to in move_to_pairs:
            for atom in atoms:
                imask += 1
                if ( not mask[imask] ):
                    continue
                distance = move_pose.residue(seqpos_move).xyz(atom).distance_squared(to_pose.residue(seqpos_to).xyz(atom))
                rmsd += distance

        rmsd /= imask
        rmsd = np.sqrt(rmsd)

        zero = numeric.xyzVector_double_t(0, 0, 0)
        xform = nup.vector_to_xform( zero - ref_com ) @ nup.matrix_to_xform( rotation_matrix ) \
                    @ nup.vector_to_xform( move_com )


    print("Final RMSD: %.3f over %i atoms"%(rmsd, mask.sum()))


    if ( return_all_move_to ):
        return rmsd, move_to_pairs, move_pose, xform, all_move_to
    else:
        return rmsd, move_to_pairs, move_pose, xform

def delete_residues_smart(pose, start, end):

    ft = pose.fold_tree()
    if ( start == 1 ):
        ft.reorder(pose.size())
    else:
        ft.reorder(1)
    pose.fold_tree(ft)

    pose.delete_residue_range_slow(start, end)
    pose.conformation().detect_disulfides()



# Uses a pareto-front style method to take the top x by multiple different values
# The general idea is to rank all values in every list (with the highest being best)
#   and then say: "If I take the top 1 best from each, do I have enough?"
#                 "If I take the top 2 best from each, do I have enough?"
#      Where when taking the top n from each, the process uses & logic, so there may actually
#          be 0 elements when taking the top 1 best from each
#
#   Each row is a design, each column is a scoreterm. Higher is better
#
#   Return value is indices to keep
def top_x_by_multiple(data, x):

    desired_num = x

    # returns indices
    tags = np.arange(len(data))

    # shuffle data so that in perfect ties, the output is a random subsample
    to_shuf = np.arange(0, len(tags), dtype=np.uint32)
    np.random.shuffle(to_shuf)
    tags = tags[to_shuf]
    data = data[to_shuf]
    total = len(data)

    # sort of like argsort for each scoreterm
    ranked = scipy.stats.rankdata(data, axis=0)

    bigger_2 = 0
    for i in range(1000):
        bigger_2 = 2**i
        if ( bigger_2 // 2 > total):
            break

    percentiles = np.linspace(0, 1, bigger_2)

    space_size = bigger_2
    next_cut = bigger_2 // 2 - 1

    remaining = 0
    cutoff = total

    last_mask_above_threshold = np.ones(total, bool)

    # binary search
    # This should never reach that many but it's better than while True
    for i in range(10000):


        # this is the actual ranking process
        # take top X in each argsort and make sure they're in all top Xs
        eval_percentile = (percentiles[next_cut] + percentiles[next_cut+1]) / 2
        cutoff = eval_percentile * total

        mask = np.ones(total, bool)

        for icol in range(data.shape[-1]):
            mask &= ranked[:,icol] >= cutoff


        remaining = mask.sum()

        if ( remaining >= desired_num ):
            last_mask_above_threshold = mask
        if ( remaining == desired_num ):
            break

        space_size //= 2

        if ( space_size == 1 ):
            break

        if ( remaining < desired_num ):
            next_cut -= space_size // 2
        else:
            next_cut += space_size // 2

    mask = last_mask_above_threshold

    keep_tags = tags[mask]
    # keep_data = data[mask]
    keep_ranked = ranked[mask]


    # not a perfect ranking. It's why we don't just use this at the start
    silly_score = keep_ranked.sum(axis=-1)

    # We need to do this because in extreme ties, we can output way more than we wanted
    final_sort = np.argsort(-silly_score)[:desired_num]

    
    return keep_tags[final_sort]


# rosetta/main/source/src/core/select/util/SelectResiduesByLayer.cc
def sidechain_neighbors(binder_Ca, binder_Cb, else_Ca):

    conevect = binder_Cb - binder_Ca
    conevect /= np.sqrt(np.sum(np.square(conevect), axis=-1))[:,None]

    vect = else_Ca[:,None] - binder_Cb[None,:]
    vect_lengths = np.sqrt(np.sum(np.square(vect), axis=-1))
    vect_normalized = vect / vect_lengths[:,:,None]

    dist_term = 1 / ( 1 + np.exp( vect_lengths - 9  ) )

    angle_term = (((conevect[None,:] * vect_normalized).sum(axis=-1) + 0.5) / 1.5).clip(0, None)

    sc_neigh = (dist_term * np.square( angle_term )).sum(axis=0)

    return sc_neigh






the_locals = None

def looping_filters(pose):

    out_score_map = {}

    pose = pose.split_by_chain()[1]
    npose = nup.npose_from_pose(pose)

    dssp = better_dssp3(pose)
    dssp = dssp.replace("E", "L")
    ss_elems = get_ss_elements2(dssp)

    longest_loop = 0
    for h, start, end in ss_elems:
        if ( h != "L" ):
            continue
        size = end - start + 1
        longest_loop = max(size, longest_loop)

    out_score_map['longest_loop'] = longest_loop


    ca = nu.extract_atoms(npose, [nu.CA])[:,:3]
    cb = nu.extract_atoms(npose, [nu.CB])[:,:3]

    sc_neigh = sidechain_neighbors( ca, cb, ca )

    # is_core = sc_neigh > 5.2
    is_core = sc_neigh > 4.9

    out_score_map['mean_scn'] = sc_neigh.mean()
    out_score_map['percent_core_scn'] = is_core.mean()




    hits, froms, tos, misses = motif_stuff2.motif_score_npose( npose )

    froms = np.array(froms)+1
    tos = np.array(tos)+1

    helix_elems = list(x for x in ss_elems if x[0] == "H")

    micro_helices = 0
    for _, start, end in helix_elems:
        size = end - start + 1
        if ( size < 9 ):
            micro_helices += 1
            # print("color red, resi %i-%i"%(start, end))
    out_score_map['micro_helices'] = micro_helices


    im_in_helix = np.zeros(pose.size()+1, int)
    im_in_helix[:] = -1
    for ihelix, (_, start, end) in enumerate(helix_elems):
        im_in_helix[start:end+1] = ihelix

    is_helix = im_in_helix > -1


    scores = []
    starts = []
    any_core_in_span = []
    window_size = 9
    for start in range(1, pose.size()+1):
        end = start + window_size - 1
        if ( end > pose.size() ):
            continue

        if ( not is_helix[start:end+1].all() ):
            continue

        our_ihelix = im_in_helix[start]
        if ( not (our_ihelix == im_in_helix[start:end+1]).all() ):
            continue


        from_us = (froms >= start) & (froms <= end)
        to_us = (tos >= start) & (tos <= end)

        interesting = from_us ^ to_us

        who = np.zeros(len(from_us), int)
        who[:] = -2
        who[~from_us] = im_in_helix[froms[~from_us]]
        who[~to_us] = im_in_helix[tos[~to_us]]

        who = who[interesting]
        who = who[who != -1]
        who = who[who != our_ihelix]
        assert( not np.any(who == -2) )

        unique_who, unique_counts = np.unique(who, return_counts=True)
        argsort = np.argsort(unique_counts)[::-1]
        sorted_who = unique_who[argsort]
        sorted_counts = unique_counts[argsort]

        hits_to_others = 0
        if ( len(sorted_counts) >= 2 ):
            hits_to_others = sorted_counts[1:].sum()


        # print(who, hits_to_others)ls

        scores.append(hits_to_others)
        starts.append(start)

        if ( is_core[start-1:end].any() ):
            any_core_in_span.append(1)
        else:
            any_core_in_span.append(0)
            # print("color yellow, resi %i-%i"%(start, end))
            # print("# %s"%(",".join("%.2f"%(x) for x in sc_neigh[start-1:end])))

    scores = np.array(scores)
    starts = np.array(starts)

    out_score_map['other_hits_9'] = (scores >= 1).mean()
    out_score_map['any_core_9'] = np.mean(any_core_in_span)




    return out_score_map





def pyro():
    return pyrosetta
def ros():
    return pyrosetta.rosetta


# don't base these on xmls so that there's no loading time
class RosettaPacker:

    def __init__(self):

        self.chainA = ros().core.select.residue_selector.ChainSelector("A")
        self.chainB = ros().core.select.residue_selector.ChainSelector("B")
        self.interface_on_A = ros().core.select.residue_selector.NeighborhoodResidueSelector(self.chainB, 10.0, False)
        self.interface_on_B = ros().core.select.residue_selector.NeighborhoodResidueSelector(self.chainA, 10.0, False)
        self.AB_interface = ros().core.select.residue_selector.OrResidueSelector( self.interface_on_A, self.interface_on_B )
        self.Not_interface = ros().core.select.residue_selector.NotResidueSelector( self.AB_interface )
        self.chainA_not_interface = ros().core.select.residue_selector.AndResidueSelector( self.Not_interface, self.chainA )
        self.chainB_not_interface = ros().core.select.residue_selector.AndResidueSelector( self.Not_interface, self.chainB )

        self.scorefxn_insta = pyro().get_fa_scorefxn()
        for term in self.scorefxn_insta.get_nonzero_weighted_scoretypes():
            name = ros().core.scoring.name_from_score_type(term)

            if ( "_dun" in name ):
                continue
            if ( "rama" in name ):
                continue
            if ( "p_aa_pp" in name ):
                continue
            if ( "fa_rep" in name ):
                continue
            if ( "fa_atr" in name ):
                continue
            if ( "fa_sol" in name ):
                continue
            if ( "hbond" in name ):
                continue
            if ( "pro_close" in name ):
                continue

            self.scorefxn_insta.set_weight(term, 0)

        self.scorefxn_insta_soft = self.scorefxn_insta.clone()
        self.scorefxn_insta_soft.set_weight(ros().core.scoring.fa_rep, 0.15)


        self.scorefxn_none = ros().core.scoring.ScoreFunctionFactory.create_score_function("none")
        self.scorefxn_atr = ros().core.scoring.ScoreFunctionFactory.create_score_function("none")
        self.scorefxn_atr.set_weight(ros().core.scoring.fa_atr, 1)
        self.scorefxn_beta = pyro().get_fa_scorefxn()
        self.scorefxn_beta_soft = ros().core.scoring.ScoreFunctionFactory.create_score_function("beta_nov16_soft")


    # pack with only dunbrack, vdw, and hbonds
    #  elec is super slow so we can't have that
    def insta_pack(self, pose):

        tf = ros().core.pack.task.TaskFactory()
        tf.push_back( ros().core.pack.task.operation.RestrictToRepacking() )
        tf.push_back( ros().core.pack.task.operation.OperateOnResidueSubset( ros().core.pack.task.operation.PreventRepackingRLT(),
                          self.chainB_not_interface, False ))
        tf.push_back( ros().core.pack.task.operation.IncludeCurrent() )

        packer = ros().protocols.minimization_packing.PackRotamersMover()
        packer.score_function( self.scorefxn_insta )
        packer.task_factory( tf )

        packer.apply( pose )

    def beta_pack(self, pose, soft=False):

        tf = ros().core.pack.task.TaskFactory()
        tf.push_back( ros().core.pack.task.operation.RestrictToRepacking() )
        tf.push_back( ros().core.pack.task.operation.OperateOnResidueSubset( ros().core.pack.task.operation.PreventRepackingRLT(),
                          self.chainB_not_interface, False ))
        tf.push_back( ros().core.pack.task.operation.IncludeCurrent() )

        packer = ros().protocols.minimization_packing.PackRotamersMover()
        if ( soft ):
            packer.score_function( self.scorefxn_beta_soft )
        else:
            packer.score_function( self.scorefxn_beta )
        packer.task_factory( tf )

        packer.apply( pose )


    # use the packer to do this to prevent issues with disulfides and to ensure we get variant types right
    def thread_seq(self, pose, new_seq):

        old_seq = pose.sequence()

        locked_subset = ros().utility.vector1_bool( pose.size() )

        tf = ros().core.pack.task.TaskFactory()
        for seqpos in range(1, pose.size()+1 ):
            old_letter = old_seq[seqpos-1]
            new_letter = new_seq[seqpos-1]

            if ( old_letter == new_letter ):
                locked_subset[seqpos] = True
                continue

            restrict_aa = ros().core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
            restrict_aa.aas_to_keep( new_letter )

            subset = ros().utility.vector1_bool( pose.size() )
            subset[seqpos] = True
            tf.push_back( ros().core.pack.task.operation.OperateOnResidueSubset( restrict_aa, subset ) )

        tf.push_back( ros().core.pack.task.operation.OperateOnResidueSubset( ros().core.pack.task.operation.PreventRepackingRLT(),
                          locked_subset ) )

        packer = ros().protocols.minimization_packing.PackRotamersMover()
        packer.score_function( self.scorefxn_none )
        packer.task_factory( tf )

        packer.apply( pose )


    def fast_hydrophobic_interface(self, pose):

        tf = ros().core.pack.task.TaskFactory()
        tf.push_back( ros().core.pack.task.operation.OperateOnResidueSubset( ros().core.pack.task.operation.PreventRepackingRLT(),
                          self.chainB, False ))
        tf.push_back( ros().core.pack.task.operation.OperateOnResidueSubset( ros().core.pack.task.operation.PreventRepackingRLT(),
                          self.chainA_not_interface, False ))
        tf.push_back( ros().core.pack.task.operation.IncludeCurrent() )

        restrict_aa = ros().core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
        restrict_aa.aas_to_keep( "GAFILMVW" )
        tf.push_back( ros().core.pack.task.operation.OperateOnResidueSubset( restrict_aa,
                          self.interface_on_A, False ))

        packer = ros().protocols.minimization_packing.PackRotamersMover()
        packer.score_function( self.scorefxn_insta_soft )
        packer.task_factory( tf )

        packer.apply( pose )


rosetta_packer = RosettaPacker()

def pose_from_silent_lines(structure, tag):
    vec = utility.vector1_std_string()
    vec.append(tag)

    stream = std.istringstream(structure)

    sfd = core.io.silent.SilentFileData(rosetta.core.io.silent.SilentFileOptions())
    sfd.read_stream(stream, vec, True, "fake")

    pose = core.pose.Pose()
    sfd.get_structure(tag).fill_pose(pose)

    return pose

def calc_ss_schain_atr(pose):


    pose = pose.split_by_chain()[1]


    size = pose.size()+1


    

    full_score = scorefxn_fa_atr(pose)

    bb_bb_scores = np.zeros((size, size))
    bb_sc_scores = np.zeros((size, size))
    sc_sc_scores = np.zeros((size, size))


    for seqpos1 in range(1, size):
        for seqpos2 in range(1, size):
            if ( seqpos1 == seqpos2 ):
                continue
            res1 = pose.residue(seqpos1)
            res2 = pose.residue(seqpos2)

            bb_bb_emap = core.scoring.EMapVector()
            scorefxn_fa_atr.eval_ci_2b_bb_bb(res1, res2, pose, bb_bb_emap)
            scorefxn_fa_atr.eval_cd_2b_bb_bb(res1, res2, pose, bb_bb_emap)

            bb_bb_scores[seqpos1, seqpos2] = bb_bb_emap.dot(scorefxn_fa_atr.weights())

            bb_sc_emap = core.scoring.EMapVector()
            scorefxn_fa_atr.eval_ci_2b_bb_sc(res1, res2, pose, bb_sc_emap)
            scorefxn_fa_atr.eval_cd_2b_bb_sc(res1, res2, pose, bb_sc_emap)

            bb_sc_scores[seqpos1, seqpos2] = bb_sc_emap.dot(scorefxn_fa_atr.weights())

            sc_sc_emap = core.scoring.EMapVector()
            scorefxn_fa_atr.eval_ci_2b_sc_sc(res1, res2, pose, sc_sc_emap)
            scorefxn_fa_atr.eval_cd_2b_sc_sc(res1, res2, pose, sc_sc_emap)

            sc_sc_scores[seqpos1, seqpos2] = sc_sc_emap.dot(scorefxn_fa_atr.weights())



    assert(np.allclose( sc_sc_scores, sc_sc_scores.T ))
    assert(np.allclose( bb_bb_scores, bb_bb_scores.T ))

    my_score = bb_bb_scores.sum()/2 + sc_sc_scores.sum()/2 + bb_sc_scores.sum()


    assert(np.isclose(my_score, full_score))


    dssp = better_dssp3(pose)
    ss_elems = get_ss_elements2(dssp)


    ss_pure_schain_atr = 0
    ss_schain_atr = 0
    ss_atr = 0
    for h1, s1, e1 in ss_elems:
        if ( h1 != "H" ):
            continue

        for h2, s2, e2 in ss_elems:
            if ( h2 != "H" ):
                continue
            if ( s1 == s2 ):
                continue

            bb_bb = bb_bb_scores[s1:e1+1,s2:e2+1].sum()/2
            bb_sc = bb_sc_scores[s1:e1+1,s2:e2+1].sum()
            sc_sc = sc_sc_scores[s1:e1+1,s2:e2+1].sum()

            ss_atr += bb_bb + bb_sc + sc_sc
            ss_schain_atr += bb_sc + sc_sc
            ss_pure_schain_atr += sc_sc

    d = {}
    d['ss_atr'] = ss_atr
    d['ss_schain_atr'] = ss_schain_atr
    d['ss_pure_schain_atr'] = ss_pure_schain_atr


    return d



def get_simple_monomer_surface(pose):
    sc_neighbors = core.select.util.SelectResiduesByLayer()
    sc_neighbors.use_sidechain_neighbors( True )
    sc_neighbors.compute(pose, "", True)

    monomer_size = pose.conformation().chain_end(1)

    is_surface = np.zeros(pose.size()+1, bool)

    for seqpos in range(1, monomer_size+1):
        sc_neigh = sc_neighbors.rsd_sasa(seqpos)

        is_surface[seqpos] = sc_neigh < 2.0

    return is_surface

def get_interface_definitely_not_monomer_core(pose):
    by_vector = interface_by_vector.apply(pose)
    small_interface = interface_on_A.apply(pose)

    pose = pose.clone()
    move_chainA_far_away(pose)

    sc_neighbors = core.select.util.SelectResiduesByLayer()
    sc_neighbors.use_sidechain_neighbors( True )
    sc_neighbors.compute(pose, "", True)

    monomer_size = pose.conformation().chain_end(1)

    is_interface = np.zeros(monomer_size+1, bool)

    for seqpos in range(1, monomer_size+1):
        sc_neigh = sc_neighbors.rsd_sasa(seqpos)

        # should be 5.2, but we're being cautious
        if sc_neigh > 4.2:
            continue
        if by_vector[seqpos] or small_interface[seqpos]:
            is_interface[seqpos] = True

    return is_interface

def calc_ddg_norepack(pose, scorefxn):
    pose = pose.clone()
    
    close_score = scorefxn(pose)
    pose = move_chainA_far_away(pose)
    far_score = scorefxn(pose)

    return close_score - far_score


def get_hbond_decoding_order(pose, interface_first=True, strict_lock=False):

    rosetta_packer.beta_pack(pose, soft=False)
    rosetta_packer.scorefxn_beta_soft(pose)

    is_interface = big_interface_on_A.apply(pose)
    is_interface_np = np.array(list(is_interface)).astype(bool)


    monomer_size = pose.conformation().chain_end(1)


    hbset = core.scoring.hbonds.HBondSet()
    core.scoring.hbonds.fill_hbond_set(pose, False, hbset)
    hbset.hbond_options().bb_donor_acceptor_check(False)
    core.scoring.hbonds.fill_hbond_set(pose, False, hbset)

    hbond_locs = []
    hb_map = np.zeros((pose.size()+1, pose.size()+1), int)

    for seqpos in range(1, monomer_size+1):

        for hbond in hbset.residue_hbonds(seqpos):

            # bb-bb hbonds don't count
            if hbond.acc_atm_is_backbone() and hbond.don_hatm_is_backbone():
                continue

            we_are_acc = hbond.acc_res() == seqpos

            we_are_bb = False
            if we_are_acc and hbond.acc_atm_is_backbone():
                we_are_bb = True

            if not we_are_acc and hbond.don_hatm_is_backbone():
                we_are_bb = True

            other_seqpos = hbond.don_res() if we_are_acc else hbond.acc_res()

            # if other_seqpos <= seqpos:
            #     continue

            if hbond.energy() > -0.2:
                continue

            don_xyz = nup.from_vector( pose.residue( hbond.don_res() ).xyz( hbond.don_hatm() ) )
            acc_xyz = nup.from_vector( pose.residue( hbond.acc_res() ).xyz( hbond.acc_atm() ) )

            hbond_locs.append( (don_xyz + acc_xyz) / 2 )

            # if we_are_bb, then we only store the location of the hbond (and not the identities)
            if we_are_bb:
                continue

            hb_map[seqpos, other_seqpos] += 1
            # hb_map[other_seqpos, seqpos] += 1


    hbond_locs = np.array(hbond_locs).reshape(-1, 3)
    hb_close_dist = 5

    close_hbs = [0]

    for seqpos in range(1, monomer_size+1):
        hb_close_mask = np.zeros(len(hbond_locs), bool)

        res = pose.residue(seqpos)

        for iatom in range(res.first_sidechain_atom(), res.nheavyatoms()+1):
            xyz = nup.from_vector( res.xyz(iatom) )

            d2 = np.sum( np.square( xyz - hbond_locs), axis=-1 )

            hb_close_mask |= d2 < hb_close_dist**2

        close_hbs.append(hb_close_mask.sum())

    close_hbs = np.array(close_hbs)


    cross_hbs = hb_map[:monomer_size+1,monomer_size+1:].sum(axis=-1)
    self_hbs = hb_map[:monomer_size+1,:monomer_size+1].sum(axis=-1)


    score = 1 * cross_hbs + 0.5 * self_hbs + 0.1 * (close_hbs - cross_hbs - self_hbs)


    n_self = 1
    if strict_lock:
        n_self = 2

    should_lock = (cross_hbs >= 2) | ((cross_hbs == 1) & (self_hbs >= n_self))

    final_scores = np.zeros(pose.size())
    final_scores[:monomer_size] = score[1:]
    final_scores[~is_interface_np] = 0
    final_scores[monomer_size:] += 10000
    if interface_first:
        final_scores[is_interface_np] += 1000
    else: 
        final_scores[~is_interface_np] += 1000

    final_scores += np.random.random(len(final_scores)) * 0.01

    order = np.argsort(-final_scores)

    locked = [x for x in np.where(should_lock)[0]]

    return order, locked


###########################################################################3
########## HIS stuff
#############################################################################


scorefxn_sc = core.scoring.ScoreFunctionFactory.create_score_function("none")
scorefxn_sc.set_weight(core.scoring.hbond_bb_sc, 1)
scorefxn_sc.set_weight(core.scoring.hbond_sc, 1)

scorefxn_bb = core.scoring.ScoreFunctionFactory.create_score_function("none")
scorefxn_bb.set_weight(core.scoring.hbond_sr_bb, 1)
scorefxn_bb.set_weight(core.scoring.hbond_lr_bb, 1)


def get_simple_is_core(pose):

    sc_neighbors = core.select.util.SelectResiduesByLayer()
    sc_neighbors.use_sidechain_neighbors( True )
    sc_neighbors.compute(pose, "", True)
    is_core = [None]
    for seqpos in range(1, pose.size()+1):
        is_core.append( sc_neighbors.rsd_sasa(seqpos) > 4)

    return is_core


def find_his_positions(pose):

    monomer_size = pose.conformation().chain_end(1)

    is_core = get_simple_is_core(pose)

    tf = core.pack.task.TaskFactory()

    repack_aa = core.pack.task.operation.RestrictToRepackingRLT()
    subset = chainB.apply(pose)
    tf.push_back( core.pack.task.operation.OperateOnResidueSubset( repack_aa, subset ) )

    restrict_aa = core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
    restrict_aa.aas_to_keep( 'H' )
    subset = chainA.apply(pose)
    tf.push_back( core.pack.task.operation.OperateOnResidueSubset( restrict_aa, subset ) )

    tf.push_back(core.pack.task.operation.InitializeFromCommandline())

    task = tf.create_task_and_apply_taskoperations( pose )

    scorefxn(pose)
    scorefxn.setup_for_packing( pose, task.repacking_residues(), task.designing_residues() )
    graph = core.pack.create_packer_graph( pose, scorefxn, task )

    rotsets = core.pack.rotamer_set.RotamerSets()
    rotsets.set_task( task )
    rotsets.initialize_pose_for_rotsets_creation( pose )
    rotsets.build_rotamers( pose, scorefxn, graph )


    complete_rotamer_sets = core.pack.rotamer_set.RotamerSets()
    position_had_rotset = utility.vector1_bool()

    hb_graph = core.pack.hbonds.hbond_graph_from_partial_rotsets(pose, rotsets, scorefxn_sc, scorefxn_bb, complete_rotamer_sets, position_had_rotset, -0.01)


    assert np.all(position_had_rotset)
    assert rotsets.nrotamers() == hb_graph.num_nodes()

    seqpos_has_his = np.zeros(pose.size()+1, dtype=bool)

    core_bb_acc = set()
    core_sc_acc = set()
    surf_bb_acc = set()
    surf_sc_acc = set()

    for ihbnode in range(1, hb_graph.num_nodes()+1):

        hbnode = hb_graph.get_node(ihbnode)
        seqpos = rotsets.res_for_rotamer(ihbnode)
        rotamer = rotsets.rotamer(ihbnode)


        if seqpos > monomer_size:
            continue

        it = hbnode.edge_list_begin( hb_graph )
        while it.valid():
            edge = _hbedge_from_lowmem(hb_graph, it.dereference(), ihbnode)
            it.pre_increment()

            we_are_first = edge.get_first_node_ind() == ihbnode

            other_node = edge.get_second_node_ind() if we_are_first else edge.get_first_node_ind()
            assert other_node != ihbnode
            other_seqpos = rotsets.res_for_rotamer(other_node)
            other_rotamer = rotsets.rotamer(other_node)

            if other_seqpos <= monomer_size:
                continue

            for hbond in edge.hbonds():

                we_are_donor = not (we_are_first ^ hbond.first_node_is_donor())


                our_iatom = hbond.local_atom_id_A()
                their_iatom = hbond.local_atom_id_D()

                if rotamer.atom_is_backbone(our_iatom):
                    continue

                their_bb = other_rotamer.atom_is_backbone(their_iatom)
                we_core = is_core[seqpos]

                if their_bb and we_core:
                    core_bb_acc.add(seqpos)
                if their_bb and not we_core:
                    surf_bb_acc.add(seqpos)
                if not their_bb and we_core:
                    core_sc_acc.add(seqpos)
                if not their_bb and not we_core:
                    surf_sc_acc.add(seqpos)


    labels_and_sets = [
    ('core_bb_acc', core_bb_acc),
    ('surf_bb_acc', surf_bb_acc),
    ('core_sc_acc', core_sc_acc),
    ('surf_sc_acc', surf_sc_acc),
    ]

    records = []

    for seqpos in range(1, monomer_size+1):
        d = {'seqpos':seqpos}
        labels = []
        for label, sett in labels_and_sets:
            d[label] = seqpos in sett
            if seqpos in sett:
                labels.append(label)

        if len(labels) > 0:
            print(f"Seqpos {seqpos:3d}: PDB {pose.pdb_info().number(seqpos):3d}: {' '.join(labels)}")

            records.append(d)

    df = pd.DataFrame(records)

    return df


def protonate_histidines(pose, do_protonate=True):

    pH = 0 if do_protonate else 14

    basic.options.set_boolean_option('pH:pH_mode', True)
    basic.options.set_real_option('pH:value_pH', pH)

    scorefxn_pH = core.scoring.ScoreFunctionFactory.create_score_function("none")
    scorefxn_pH.set_weight(core.scoring.e_pH, 100)

    his_sel = core.select.residue_selector.ResidueNameSelector()
    his_sel.set_residue_name3("HIS")

    his_sub = his_sel.apply(pose)

    protocols.toolbox.pose_manipulation.repack_these_residues(his_sub, pose, scorefxn_pH)

    protonated_his_check(pose, do_protonate)

    basic.options.set_boolean_option('pH:pH_mode', False)
    basic.options.set_real_option('pH:value_pH', 7)


def protonated_his_check(pose, do_protonate):
    his_sel = core.select.residue_selector.ResidueNameSelector()
    his_sel.set_residue_name3("HIS")

    his_sub = his_sel.apply(pose)
    for seqpos in range(1, pose.size()+1):
        if not his_sub[seqpos]:
            continue
        if do_protonate:
            assert '_P' in pose.residue(seqpos).name(), 'Protonation failed. Do you have -pH_mode True?'
        else:
            assert '_P' not in pose.residue(seqpos).name(), 'Deprotonation failed. Do you have -pH_mode True?'


def drop_reslabel(pose, label_re):
    compiled = re.compile(label_re)
    for seqpos in range(1, pose.size()+1):
        labels = pose.pdb_info().get_reslabels(seqpos)
        if ( len(labels) == 0 ):
            continue
        keep_mask = []
        for label in labels:
            if ( compiled.match(label) ):
                keep_mask.append(False)
            else:
                keep_mask.append(True)
        if ( not np.all(keep_mask) ):
            pose.pdb_info().clear_reslabel(seqpos)
            for label, keep in zip(labels, keep_mask):
                if ( not keep ):
                    continue
                pose.pdb_info().add_reslabel(seqpos, label)


scorefxn_hbonds = get_fa_scorefxn()
fix_scorefxn(scorefxn_hbonds, allow_double_bb=True)


def find_good_his_ph_interactions(in_pose):
    pose = in_pose.clone()

    monomer_size = pose.conformation().chain_end(1)

    is_core = get_simple_is_core(pose)


    big_df = None

    for low_pH in [False, True]:

        records = []

        protonate_histidines(pose, low_pH)
        rosetta_packer.beta_pack(pose, soft=True)

        protonated_his_check(pose, low_pH)

        scorefxn_hbonds(pose)
        hbset = core.scoring.hbonds.HBondSet()
        core.scoring.hbonds.fill_hbond_set(pose, False, hbset)
        hbset.hbond_options().bb_donor_acceptor_check(False)
        core.scoring.hbonds.fill_hbond_set(pose, False, hbset)

        prefix = 'low_pH_'if low_pH else 'high_pH_'

        cats = ['core_A', 'core_D', 'core_AD', 'surf_A', 'surf_D', 'surf_AD']
        if low_pH:
            cats = ['core_D', 'core_DD', 'surf_D', 'surf_DD']

        cats2 = list(cats)
        for cat in cats:
            cats2.append(cat + '_bb')


        for seqpos in range(0, pose.size()+1):
            d = {'seqpos':str(seqpos)}
            for cat in cats2:
                d[prefix + cat] = 0

            if seqpos == 0: # make sure there's a df
                records.append(d)
                continue

            if pose.residue(seqpos).name1() != 'H':
                continue

            we_are_binder = seqpos <= monomer_size

            his_ACC = 0
            his_DON = 0
            cross_hbonds = 0
            other_is_bb = False
            atom_used = set()
            for hbond in hbset.residue_hbonds(seqpos):
                we_are_don = hbond.don_res() == seqpos

                to_bb = False
                if we_are_don:
                    if hbond.don_hatm_is_backbone():
                        continue
                    atom = hbond.don_hatm()
                    to_bb = hbond.acc_atm_is_backbone()
                else:
                    if hbond.acc_atm_is_backbone():
                        continue
                    atom = hbond.acc_atm()
                    to_bb = hbond.don_hatm_is_backbone()

                # So ok, theoretically there could be a bug if the same atom is making a hbond to both the target
                #  and to the binder and the binder atom comes up first
                # But like... What are the odds?
                if atom in atom_used:
                    continue
                atom_used.add(atom)

                if we_are_don:
                    his_DON += 1
                else:
                    his_ACC += 1

                if we_are_don:
                    other_res = hbond.acc_res()
                else:
                    other_res = hbond.don_res()

                other_is_binder = other_res <= monomer_size

                if we_are_binder != other_is_binder:
                    cross_hbonds += 1
                    other_is_bb = other_is_bb or to_bb

            if cross_hbonds == 0:
                continue
            if his_ACC + his_DON == 0:
                continue


            assert his_ACC + his_DON <= 2

            # this is actually sidechain neighbors
            core_str = 'core_' if is_core[seqpos] else 'surf_'

            ad_str = "A"*his_ACC + "D"*his_DON

            cat = core_str + ad_str
            if other_is_bb:
                cat += '_bb'
            assert cat in cats2, f'How did we get {cat} with low_pH: {low_pH}'

            assert prefix + cat in d
            d[prefix + cat] = cross_hbonds
            records.append(d)

            print(f'{prefix} {seqpos:3d} {pose.pdb_info().chain(seqpos)} {cat:7s} cross hbonds: {cross_hbonds}')


        df = pd.DataFrame(records)
        if big_df is None:
            big_df = df
        else:
            print(list(big_df))
            big_df = df.merge(big_df, 'outer', 'seqpos')

    big_df['seqpos'] = big_df['seqpos'].astype(int)
    big_df = big_df[big_df['seqpos'] > 0]
    big_df = big_df.fillna(0)

    big_df['ph_score'] = calc_ph_score(big_df)

    out_pose = in_pose.clone()
    drop_reslabel(out_pose, 'ph_score.*')
    for idx, row in big_df.iterrows():
        score = row['ph_score']
        out_pose.pdb_info().add_reslabel(int(row['seqpos']), 'ph_score:%.1f'%score)

    # print(big_df[['ph_score']])

    return big_df, out_pose





def calc_ph_score(df):

    return (
       -1.0 * df['low_pH_surf_D'] + 
       -3.0 * df['low_pH_surf_DD'] + 
       -2.0 * df['low_pH_surf_D_bb'] + 
       -6.0 * df['low_pH_surf_DD_bb'] + 
       -3.0 * df['low_pH_core_D'] + 
       -9.0 * df['low_pH_core_DD'] + 
       -6.0 * df['low_pH_core_D_bb'] + 
       -18.0 * df['low_pH_core_DD_bb'] + 

        1.0 * df['high_pH_surf_A'] + 
        1.0 * df['high_pH_surf_D'] + 
        3.0 * df['high_pH_surf_AD'] + 
        2.0 * df['high_pH_surf_A_bb'] + 
        2.0 * df['high_pH_surf_D_bb'] + 
        6.0 * df['high_pH_surf_AD_bb'] +

        3.0 * df['high_pH_core_A'] + 
        3.0 * df['high_pH_core_D'] + 
        9.0 * df['high_pH_core_AD'] + 
        6.0 * df['high_pH_core_A_bb'] + 
        6.0 * df['high_pH_core_D_bb'] + 
        18.0 * df['high_pH_core_AD_bb'] 

        )




def draw_sets(main_probs, seqposs, n_sets, groups_of, overtry=100, used_decay=0.5):

    main_probs = main_probs.copy()

    sets_to_test = set()
    for i in range(args.first_round_n_sets * overtry):
        if len(sets_to_test) == args.first_round_n_sets:
            continue
        
        probs = main_probs.copy()
        probs /= probs.sum()

        n_draw = min(len(main_probs), args.his_in_groups_of)

        my_set = tuple(list(sorted(list(np.random.choice(seqposs, size=n_draw, replace=False, p=probs)))))

        if my_set in sets_to_test:
            continue

        sets_to_test.add(my_set)

        for seqpos in my_set:
            idx = np.where(seqposs == seqpos)[0]
            main_probs[idx] *= used_decay

    return list(sets_to_test)

def do_his_trials(pose, name_no_suffix, sets, seqpos_df):
    global force_seq_by_tag
    for key in list(force_seq_by_tag):
        del force_seq_by_tag[key]

    poses = []
    dfs = []
    position_scores = seqpos_df.copy()
    position_scores['seqpos'] = position_scores['seqpos'].astype(str)
    position_scores['best_his'] = -1

    monomer_size = pose.conformation().chain_end(1)

    for sett in sets:
        sett_str = '_'.join(['%iH'%x for x in sett])
        this_name = name_no_suffix + '_' + sett_str

        this_seq = ['.']*monomer_size
        for seqpos in sett:
            this_seq[seqpos-1] = 'H'
        force_seq_by_tag[this_name] = this_seq

        this_poses, this_score_df = get_poses_and_scores(pose, this_name, n_per_input=1)
        if len(this_poses) == 0:
            continue
        assert len(this_poses) == 1

        poses.append(this_poses[0])
        dfs.append(this_score_df)

        his_df, labeled_pose = find_good_his_ph_interactions(this_poses[0][0])
        this_poses[0][0] = labeled_pose

        this_score_df['ph_score'] = his_df['ph_score'].sum()

        his_df = his_df[['ph_score', 'seqpos']].copy()
        his_df['seqpos'] = his_df['seqpos'].astype(str)
        position_scores = position_scores.merge(his_df, 'left', 'seqpos')
        position_scores['ph_score'] = position_scores['ph_score'].fillna(0)

        # The best_his gets the best his score we've ever seen
        position_scores['best_his'] = np.maximum(position_scores['best_his'], position_scores['ph_score'])
        position_scores = position_scores.drop('ph_score', axis=1)



        print("pH SCORE: %3i -- %s"%(this_score_df['ph_score'].iloc[0], sett_str))


    if len(poses) == 0:
        return [], None, position_scores

    return poses, pd.concat(dfs), position_scores


def do_ph_design(pose, name_no_suffix):
    global force_seq_by_tag


    choice_df = find_his_positions(pose)

    if len(choice_df) == 0:
        return [], None

    choice_df['initial_his_score'] = (
        6.0 * choice_df['core_bb_acc'] + 
        3.0 * choice_df['core_sc_acc'] + 
        2.0 * choice_df['surf_bb_acc'] + 
        1.0 * choice_df['surf_sc_acc']
        )


    sets = draw_sets(
        choice_df['initial_his_score'].values, 
        choice_df['seqpos'].values,
        args.first_round_n_sets,
        args.his_in_groups_of,
        )

    print('First round testing:')
    print('\n'.join([' '.join('%3i'%i for i in x) for x in sets]))
    first_poses, first_df, position_scores = do_his_trials(pose, name_no_suffix, sets, choice_df[['seqpos']])

    if len(first_poses) == 0:
        print("There were no outputs?")
        return [], None

    choice_df['seqpos'] = choice_df['seqpos'].astype(str)
    choice_df = choice_df.merge(position_scores, 'left', 'seqpos')
    choice_df['seqpos'] = choice_df['seqpos'].astype(int)

    choice_df = choice_df[choice_df['best_his'] > 0]

    if len(choice_df) == 0:
        print("None of the histidines were any good")
        return [], None

    print('Histidine scores')
    for idx, row in choice_df.sort_values('best_his', ascending=False).iterrows():
        print('%2iH: %3i'%(row['seqpos'], row['best_his']))

    sets = draw_sets(
        choice_df['best_his'].values, 
        choice_df['seqpos'].values,
        args.second_round_n_sets,
        args.his_in_groups_of,
        )

    print('Second round testing:')
    print('\n'.join([' '.join('%3i'%i for i in x) for x in sets]))
    second_poses, second_df, _ = do_his_trials(pose, name_no_suffix, sets, choice_df[['seqpos']])

    if len(second_poses) > 0:
        first_poses += second_poses
        first_df = pd.concat((first_df, second_df))
        first_df = first_df.reset_index(drop=True)


    final_mask = first_df['ph_score'] > args.min_his_score
    if final_mask.sum() > args.num_per_input:
        wh_final_mask = np.where(final_mask)[0]
        those_scores = first_df['ph_score'].values[final_mask]
        that_argsort = np.argsort(-those_scores)
        keep_args = that_argsort[:args.num_per_input]
        wh_keep = wh_final_mask[keep_args]
        final_mask[:] = False
        final_mask[wh_keep] = True

    print("Found %i outputs"%final_mask.sum())

    if final_mask.sum() == 0:
        return [], None

    final_poses = [first_poses[x] for x in np.where(final_mask)[0]]
    final_df = first_df[final_mask]

    return final_poses, final_df

######################################################################



the_locals = None

def get_poses_and_scores(in_pose, name_no_suffix, n_per_input=1):

    if ( args.do_monomer_filter ):
        looping_results = looping_filters(in_pose)

        if ( looping_results['longest_loop'] > args.longest_loop ):
            print("Fail longest_loop", looping_results['longest_loop'])
            return [],[]
        if ( looping_results['percent_core_scn'] < args.percent_core_scn ):
            print("Fail percent_core_scn", looping_results['percent_core_scn'])
            return [],[]
        if ( looping_results['micro_helices'] > args.micro_helices ):
            print("Fail micro_helices", looping_results['micro_helices'])
            return [],[]
        if ( looping_results['any_core_9'] < args.any_core_9 ):
            print("Fail any_core_9", looping_results['any_core_9'])
            return [],[]
        if ( looping_results['other_hits_9'] < args.other_hits_9 ):
            print("Fail other_hits_9", looping_results['other_hits_9'])
            return [],[]

    monomer_size = in_pose.conformation().chain_end(1)

    # delete stuff that's far away from the terminis to make the rest faster
    if ( args.crappy_interface_mode ):
        npose = nup.npose_from_pose(in_pose)

        cas = nu.extract_atoms(npose, [nu.CA])[:,:3]

        realized_patch_idx = np.array(list(all_patch_indices)) + monomer_size - 1

        patch_cas = cas[realized_patch_idx]
        monomer_cas = cas[:monomer_size]

        monomer_to_patch = np.linalg.norm( monomer_cas[:,None] - patch_cas[None,:], axis=-1 )
        closest_to_patch = monomer_to_patch.min(axis=-1)

        close_enough = closest_to_patch < args.trim_distance
        if ( not close_enough.any()):
            return [],[]

        last_removed_res_front_idx1 = list(close_enough).index(True)
        first_removed_res_back_idx1 = monomer_size+1 - list(close_enough[::-1]).index(True)

        if ( first_removed_res_back_idx1 <= monomer_size ):
            print("Removing %i residues from cterm"%(monomer_size-first_removed_res_back_idx1+1))
            delete_residues_smart(in_pose, first_removed_res_back_idx1, monomer_size)
        if ( last_removed_res_front_idx1 >= 1 ):
            print("Removing %i residues from nterm"%(last_removed_res_front_idx1-1+1))
            delete_residues_smart(in_pose, 1, last_removed_res_front_idx1)

    if ( new_target is not None ):

        my_target = new_target.clone()

        rmsd, move_to_pairs, move_pose, xform = pymol_align( my_target, in_pose, sel_move=None, sel_to=chainB.apply(in_pose),
                                to_pose_upweight_mask=interface_on_B.apply(in_pose), return_all_move_to=False )


        in_pose = in_pose.split_by_chain()[1]
        in_pose.append_pose_by_jump(my_target, 1)
        in_pose.pdb_info(core.pose.PDBInfo(in_pose))

        # in_pose.dump_pdb("test3.pdb")
        # my_target.dump_pdb("test2.pdb")
        # sys.exit()

    tied_translate = numeric.xyzVector_double_t(1000, 0, 0)
    tied_untranslate = numeric.xyzVector_double_t(-1000, 0, 0)

    tied_pose = None
    tied_positions_list = None
    if tied_target is not None:

        my_target = tied_target.clone()
        monomer_size = in_pose.conformation().chain_end(1)

        rmsd, move_to_pairs, move_pose, xform = pymol_align( my_target, in_pose, sel_move=None, sel_to=chainB.apply(in_pose),
                                to_pose_upweight_mask=interface_on_B.apply(in_pose), return_all_move_to=False )


        tied_pose = in_pose.split_by_chain()[1]
        tied_pose.append_pose_by_jump(my_target, 1)
        pdb_info = core.pose.PDBInfo(tied_pose)

        for seqpos in range(1, monomer_size+1):
            pdb_info.chain(seqpos, 'C')

        for seqpos in range(monomer_size+1, tied_pose.size()+1):
            pdb_info.chain(seqpos, 'D')


        tied_pose.pdb_info(pdb_info)

        assert in_pose.pdb_info().number(1) == 1
        assert in_pose.pdb_info().number(monomer_size) == monomer_size


        tied_positions_list = []
        for seqpos in range(1, monomer_size+1):
            tied_positions_list.append({'A':[seqpos], 'C':[seqpos]})

        tied_pose.apply_transform_Rx_plus_v(numeric.xyzMatrix_double_t.identity(), tied_translate)



    monomer_size = in_pose.conformation().chain_end(1)

############ mpnn surface bias ###################
    is_surface = get_simple_monomer_surface(in_pose)
    is_loop = np.array([x == "L" for x in better_dssp3(in_pose)])

    # is_surface = is_surface & ~is_loop

    alphabet = 'ACDEFGHIKLMNPQRSTVWYX'

    per_pos_bias = np.zeros([monomer_size, 21])
    for seqpos in range(1, monomer_size+1):
        if ( not is_surface[seqpos] or is_loop[seqpos] ):
            continue
        for il, letter in enumerate(alphabet):
            if ( letter in surface_bias_AA_dict ):
                per_pos_bias[seqpos-1,il] += surface_bias_AA_dict[letter]
                

############# interface desapping ######################
    is_interface = get_interface_definitely_not_monomer_core(in_pose)
    print('+'.join(str(x) for x in np.where(is_interface)[0]))

    for seqpos in range(1, monomer_size+1):
        if not is_interface[seqpos]:
            continue
        for il, letter in enumerate(alphabet):
            if letter in interface_bias_AA_dict:
                per_pos_bias[seqpos-1,il] += interface_bias_AA_dict[letter]


############ pssm ######################################
    aa_to_idx = {}
    for il, letter in enumerate(alphabet):
        aa_to_idx[letter] = il

    if pssms is not None:
        assert( name_no_suffix in pssms)
        pssm = eval(pssms[name_no_suffix])
        print("Found pssm")
        for seqpos in pssm:
            for letter in pssm[seqpos]:
                value = float(pssm[seqpos][letter])
                idx = aa_to_idx[letter]
                per_pos_bias[int(seqpos)-1,idx] += value * args.pssm_multiplier



    bias_by_res_d = {}
    bias_by_res_d['A'] = per_pos_bias


##################################################





    mpnn_fixed = None
    if ( len(force_seq_by_tag) != 0 ):
        assert(name_no_suffix in force_seq_by_tag)

        force_seq = force_seq_by_tag[name_no_suffix]
        assert(len(force_seq) == monomer_size)

        fixed_positions = []
        new_monomer_seq = ""
        for i in range(len(force_seq)):
            if ( force_seq[i] == "." ):
                new_monomer_seq += "A"
            else:
                new_monomer_seq += force_seq[i]
                fixed_positions.append(i+1)


        new_seq = new_monomer_seq + in_pose.sequence()[monomer_size:]

        rosetta_packer.thread_seq(in_pose, new_seq)

        mpnn_fixed = {"A":fixed_positions}



##################################################

    # hbnet locking

    if args.hbnet_lock_identities or args.hbnet_lock_identities_strict:

        this_decoding_order, this_fixed_pos = get_hbond_decoding_order(in_pose, strict_lock=args.hbnet_lock_identities_strict)

        print("HBNet Locking:", this_fixed_pos)
        if mpnn_fixed is None:
            mpnn_fixed = {"A":this_fixed_pos}
        else:
            mpnn_fixed = {"A":list(set(this_fixed_pos) | set(mpnn_fixed['A']))}





#################################################

    in_sequence = in_pose.sequence()
    temps=[float(x) for x in args.mpnn_temps.split()]


    mpnn_pose = in_pose.clone()

    if tied_pose is not None:
        mpnn_pose.append_pose_by_jump(tied_pose, 1)

        if bias_by_res_d is not None:
            bias_by_res_d['C'] = bias_by_res_d['A']

        if mpnn_fixed is not None:
            mpnn_fixed['C'] = mpnn_fixed['A']


    ss = ros().std.stringstream()
    mpnn_pose.dump_pdb(ss)
    pdb_str = ss.str()

    mpnn_seqs, mpnn_scores = mpnn(pdb_str, best_of_n=args.mpnn_seqs, return_all=True, sampling_temps=temps, bias_by_res_d=bias_by_res_d,
                                                                                                fixed_pos_list_by_chain_d=mpnn_fixed,
                                                                                                tied_positions_list=tied_positions_list
                                                                                                # N_solve_cycles=args.mpnn_n_cycles,
                                                                                                # single_chain=args.mpnn_single_chain
                                                                                                )


    print(mpnn_seqs[0])
    all_scores = []
    all_poses = []

    for iseq, mpnn_seq in enumerate(mpnn_seqs):

        if tied_pose is None:
            tied_iters = 1
            assert len(mpnn_seq) == monomer_size
        else:
            tied_iters = 2
            assert len(mpnn_seq) == 2*monomer_size
            assert mpnn_seq[:monomer_size] == mpnn_seq[monomer_size:]
            mpnn_seq = mpnn_seq[:monomer_size]

        these_scores = {}
        these_poses = []

        for itied in range(tied_iters):

            if itied == 0:
                pose = in_pose.clone()
                s = ''
            else:
                pose = tied_pose.clone()
                pose.apply_transform_Rx_plus_v(numeric.xyzMatrix_double_t.identity(), tied_untranslate)
                pose.pdb_info(core.pose.PDBInfo(pose))
                s = '_tied'

            target_sequence = pose.sequence()[monomer_size:]

            new_seq = mpnn_seq + target_sequence
            assert(len(new_seq) == pose.size())

            rosetta_packer.thread_seq(pose, new_seq)

            assert(pose.sequence() == new_seq)
            rosetta_packer.beta_pack(pose, soft=False)


            monomer, target = pose.split_by_chain()

            these_scores['sap_score'+s] = core.pack.guidance_scoreterms.sap.calculate_sap(pose, chainA, chainA, chainA)

            if ( these_scores['sap_score'+s] > args.sap_limit ):
                break

            these_scores['ss_schain_atr'+s] = calc_ss_schain_atr(monomer)['ss_schain_atr']
            these_scores['ddg_norepack_soft'+s] = calc_ddg_norepack(pose, rosetta_packer.scorefxn_beta_soft)
            these_scores['contact_patch'+s] = calc_cp_score(pose)
            print(these_scores)

            these_poses.append(pose)

        if ( these_scores['sap_score'+s] > args.sap_limit ):
            continue

        these_scores['parent'] = name_no_suffix
        these_scores['description'] = name_no_suffix + "_mpnn%03i"%(iseq+1)


        all_scores.append(these_scores)
        all_poses.append(these_poses)

    if ( len(all_scores) == 0 ):
        return [], []


    scores_df = pd.DataFrame(all_scores)

    # to_filter = np.zeros((len(scores_df), 4))

    filter_terms = ['contact_patch', 'ddg_norepack_soft', 'ss_schain_atr', 'sap_score']

    if not args.filter_on_sap_too:
        filter_terms.remove('sap_score')

    if args.crappy_interface_mode:
        filter_terms = ['contact_patch']

    if args.dont_filter_patch:
        filter_terms.remove('contact_patch')


    to_filter = []
    for itied in range(tied_iters):
        s = '' if itied == 0 else '_tied'

        for term in filter_terms:
            if term == 'contact_patch' and itied == 0: # contact patch res are probably wrong on the other side
                to_filter.append( scores_df['contact_patch'+s] )
            if term == 'ddg_norepack_soft':
                to_filter.append( -scores_df['ddg_norepack_soft'+s] )
            if term == 'ss_schain_atr':
                to_filter.append( -scores_df['ss_schain_atr'+s] )
            if term == 'sap_score':
                to_filter.append( -(scores_df['sap_score'+s].clip(25, None)))

    to_filter = np.stack(to_filter, axis=-1)

    to_filter = np.zeros((len(scores_df), 4))

    # to_filter[:,0] = scores_df['contact_patch']
    # to_filter[:,1] = -scores_df['ddg_norepack_soft']
    # to_filter[:,2] = -scores_df['ss_schain_atr']
    # to_filter[:,3] = -(scores_df['sap_score'].clip(25, None))

    # if ( not args.filter_on_sap_too ):
    #     to_filter = to_filter[:,:3]


    # # we only contact patch in crappy interface mode
    # if ( args.crappy_interface_mode ):
    #     to_filter = to_filter[:,:1]

    # if ( args.dont_filter_patch ):
    #     to_filter = to_filter[:,1:]

    # assert np.allclose(to_filter2, to_filter)

    keep_indices = top_x_by_multiple(to_filter, n_per_input)

    final_df = scores_df.iloc[keep_indices].copy()
    final_poses = [all_poses[x] for x in keep_indices]


    return final_poses, final_df 






############### BEGIN MAIN FUNCTION ###########################

if ( silent != "" ):
    sfd_in = rosetta.core.io.silent.SilentFileData(rosetta.core.io.silent.SilentFileOptions())
    sfd_in.read_file(silent)

    pdbs = list(sfd_in.tags())

    sfd_out = core.io.silent.SilentFileData( "out.silent", False, False, "binary", core.io.silent.SilentFileOptions())



ckpt = 0
if (os.path.exists("ckpt")):
    with open("ckpt") as f:
        try:
            ckpt = int(f.read().strip())
        except:
            pass




all_poses = []
all_scores = []

num = -1
for ipdb, pdb in enumerate(pdbs):
    t0 = time.time()
    print("Attempting pose: " + pdb)


    if ( ckpt > ipdb ):
        print("Checkpoint: Already finished: " + pdb)
        continue
    with open("ckpt", "w") as f:
        f.write(str(ipdb))

    # try:
    for k in [1]:
        if ( silent == "" ):
            pose = pose_from_file(pdb)
        else:
            pose = Pose()
            sfd_in.get_structure(pdb).fill_pose(pose)

        name_no_suffix = my_rstrip(my_rstrip(os.path.basename(pdb), ".gz"), ".pdb")

        sfd = core.io.raw_data.ScoreFileData("score.sc")

        score_map = std.map_std_string_double()
        string_map = std.map_std_string_std_string()


        final_poses, final_df = do_ph_design(pose, name_no_suffix)

        for i_final in range(len(final_poses)):
            poses = final_poses[i_final]
            for itied in range(len(poses)):
                pose = poses[itied]
                scores = final_df.iloc[i_final]

                score_map = std.map_std_string_double()
                string_map = std.map_std_string_std_string()

                description = scores['description'] if itied == 0 else scores['description'] + '_tied'

                for score in list(final_df):
                    if ( score == "description" ):
                        continue
                    if ( score == "parent" ):
                        continue

                    score_map[score] = scores[score]

                sfd = core.io.raw_data.ScoreFileData("score.sc")
                sfd.write_pose( pose, score_map, description, string_map)

                silent_name = "out.silent"
                sfd_out = core.io.silent.SilentFileData( silent_name, False, False, "binary", core.io.silent.SilentFileOptions())
                struct = sfd_out.create_SilentStructOP()
                struct.fill_struct(pose, description)
                for score in score_map:
                    struct.add_energy(score, score_map[score])
                sfd_out.add_structure(struct)
                sfd_out.write_all(silent_name, False)



        del final_poses
        del final_df


        seconds = int(time.time() - t0)

        print("protocols.jd2.JobDistributor: " + name_no_suffix + " reported success in %i seconds"%seconds)

        # except Exception as e:




with open("ckpt", "w") as f:
    f.write(str(len(pdbs)))















