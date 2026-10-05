#!/usr/bin/env python
from __future__ import division

# Exploder step 1: finds h-bonding HIS pairs (or HIS next to R/K) that can be placed in the binder.
#
# Usage: ./his_ph_exploder_v3.py pdb1.pdb pdb2.pdb [--other_aas RK] [--interface_mode] [--two_sided_design] [--dist_cutoff_mode A] [--only_allow_positions 1,2,3]
#    or: ./his_ph_exploder_v3.py -in:file:silent my.silent [--other_aas RK] [--interface_mode] [--two_sided_design] [--dist_cutoff_mode A] [--only_allow_positions 1,2,3]

import os
import sys
import math

import os
import sys
#import silent_tools

from pyrosetta import *
from pyrosetta.rosetta import *

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

import pandas as pd

def _hbedge_from_lowmem(hb_graph, lowmem_edge, node_ind):
    # The edge list iterates over LowMemEdges. We need the full HBondEdge (with the hbonds).
    if hasattr(hb_graph, 'HBondEdge_from_LowMemEdge'):  # pyrosetta with the hbond graph patch
        return hb_graph.HBondEdge_from_LowMemEdge(lowmem_edge)
    # Published pyrosetta: look the full edge up from its two node indices
    return hb_graph.find_edge(node_ind, lowmem_edge.get_other_ind(node_ind))

# import pyRMSD.RMSDCalculator

init("-beta_nov16 -in:file:silent_struct_type binary -keep_input_scores false -mute all -ex1 -ex2"
    " -holes:dalphaball /work/tlinsky/Rosetta/main/source/external/DAlpahBall/DAlphaBall.macgcc"
    )




parser = argparse.ArgumentParser()
parser.add_argument("-in:file:silent", type=str, default="")
parser.add_argument("pdbs", type=str, nargs="*")
parser.add_argument("--other_aas", default='')
parser.add_argument("--only_allow_positions", default='', help='comma separated list of positions 1-indexed we can mutate')
parser.add_argument("--interface_mode", action='store_true', help='Only allow h-bonds that cross the A/B interface')
parser.add_argument("--two_sided_design", action='store_true', help='Allow the target to mutate')
parser.add_argument("--dist_cutoff_mode", type=float, default=-1, help='Instead of looking for h-bonds, '
                                                                            'the sidechains must have non-C atoms within this many A')


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
big_interface_on_B = core.select.residue_selector.NeighborhoodResidueSelector(chainA, 14.0, False)
interface_by_vector = core.select.residue_selector.InterGroupInterfaceByVectorSelector(interface_on_A, interface_on_B)
interface_by_vector.cb_dist_cut(11)
interface_by_vector.cb_dist_cut(5.5)
interface_by_vector.vector_angle_cut(75)
interface_by_vector.vector_dist_cut(9)


A_or_big_B = core.select.residue_selector.OrResidueSelector(chainA, big_interface_on_B)



    
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
        to_pose_upweight_mask = np.zeros(to_pose.size()+1, np.bool)

    if ( move_pose_upweight_mask is None ):
        move_pose_upweight_mask = np.zeros(move_pose.size()+1, np.bool)

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


def pose_from_silent_lines(structure, tag):
    vec = utility.vector1_std_string()
    vec.append(tag)

    stream = std.istringstream(structure)

    sfd = core.io.silent.SilentFileData(rosetta.core.io.silent.SilentFileOptions())
    sfd.read_stream(stream, vec, True, "fake")

    pose = core.pose.Pose()
    sfd.get_structure(tag).fill_pose(pose)

    return pose



def classify_binder_positions(pose, interface_mode=False):

    monomer = pose.split_by_chain()[1]
    sequence = monomer.sequence()

    dssp = better_dssp3(monomer)

    if interface_mode:
        target = pose.split_by_chain()[2]
        sequence += target.sequence()
        dssp += better_dssp3(target)[1:]

    sc_neighbors = core.select.util.SelectResiduesByLayer()
    sc_neighbors.use_sidechain_neighbors( True )
    sc_neighbors.compute(pose, "", True)

    atomic_depth = core.scoring.atomic_depth.AtomicDepth( pose, 2.3, False, 0.5 )
    atomic_depth_monomer = core.scoring.atomic_depth.AtomicDepth( monomer, 2.3, False, 0.5 )
    if interface_mode:
        atomic_depth_target = core.scoring.atomic_depth.AtomicDepth( target, 2.3, False, 0.5 )
    type_set = pose.residue(1).type().atom_type_set()

    probe_size = 2.8
    per_atom_sasa = core.id.AtomID_Map_double_t()
    rsd_sasa = utility.vector1_double()
    core.scoring.calc_per_atom_sasa(pose, per_atom_sasa, rsd_sasa, 2.8, False)

    scorefxn(pose)
    scorefxn(monomer)
    if interface_mode:
        scorefxn(target)

    interface_subset = interface_by_vector.apply(pose)

    pose_dats = []

    for seqpos in range(1, pose.size()+1 if interface_mode else monomer.size()+1):

        data = {"description":name_no_suffix, "seqpos":seqpos}

        is_target = seqpos > monomer.size()

        data['sc_neighbors'] = sc_neighbors.rsd_sasa(seqpos)
        data['is_loop'] = dssp[seqpos] == "L"
        data['by_vector'] = interface_subset[seqpos]
        data['dssp'] = dssp[seqpos]
        data['name1'] = pose.residue(seqpos).name1()

        res = pose.residue(seqpos)
        if is_target:
            monomer_res = target.residue(seqpos-monomer.size())
        else:
            monomer_res = monomer.residue(seqpos)
        data['depth'] = atomic_depth.calcdepth(res.atom(res.nbr_atom()), type_set)
        if is_target:
            data['depth_monomer'] = atomic_depth_target.calcdepth(monomer_res.atom(monomer_res.nbr_atom()), type_set)
        else:
            data['depth_monomer'] = atomic_depth_monomer.calcdepth(monomer_res.atom(monomer_res.nbr_atom()), type_set)

        if is_target:
            data['ddg'] = 2*(pose.energies().residue_total_energy(seqpos) - target.energies().residue_total_energy(seqpos-monomer.size()))
        else:
            data['ddg'] = 2*(pose.energies().residue_total_energy(seqpos) - monomer.energies().residue_total_energy(seqpos))

        sc_sasa = 0
        for i in range(res.first_sidechain_atom(), res.nheavyatoms()+1):
            sc_sasa += per_atom_sasa(seqpos, i)
            for j in range(res.attached_H_begin(i), res.attached_H_end(i)+1):
                sc_sasa += per_atom_sasa(seqpos, j)

        data['sc_sasa'] = sc_sasa

        pose_dats.append(data)

    # sequences[name] = sequence

    pose_df = pd.DataFrame(pose_dats)

    # These should be option flags. But this calculation is so complicated that unless you're looking at the code, you're not going
    #  to specify this stuff correctly.
    # So either change the hardcoded stuff, or make the option flags yourself. (Or even change the calculation here)
    has_ddg = np.abs(pose_df['ddg']) > 1
    is_core = ((pose_df['sc_neighbors'] > 5.2) | (pose_df['depth'] - pose_df['depth_monomer'] > 1 )) & (pose_df['sc_sasa'] < 5)

    pose_df['is_interface_core'] = has_ddg & is_core
    pose_df['is_interface_boundary'] = ( has_ddg | pose_df['by_vector'] ) & ~pose_df['is_interface_core']
    pose_df['is_monomer_core'] = ~has_ddg & (pose_df['sc_neighbors'] >= 5.2) & ~pose_df['is_interface_boundary']
    pose_df['is_monomer_boundary'] = ~has_ddg & (pose_df['sc_neighbors'] < 5.2) & (pose_df['sc_neighbors'] >= 2.0) & ~pose_df['by_vector']
    pose_df['is_monomer_surface'] = ~has_ddg & (pose_df['sc_neighbors'] < 2.0) & ~pose_df['by_vector']

    # Make sure all positions are in exactly 1 category
    assert( np.all( pose_df[['is_interface_core', 'is_interface_boundary', 'is_monomer_core', 'is_monomer_boundary',
                            'is_monomer_surface']].sum(axis=1) == 1) )

    pose_df['seqpos'] = pose_df['seqpos'].astype(str)
    pose_df['seqpos_int'] = pose_df['seqpos'].astype(int)


    pose_df['class_label'] = 'None'
    pose_df.loc[pose_df['is_interface_core'], 'class_label'] = 'is_interface_core'
    pose_df.loc[pose_df['is_interface_boundary'], 'class_label'] = 'is_interface_boundary'
    pose_df.loc[pose_df['is_monomer_core'], 'class_label'] = 'is_monomer_core'
    pose_df.loc[pose_df['is_monomer_boundary'], 'class_label'] = 'is_monomer_boundary'
    pose_df.loc[pose_df['is_monomer_surface'], 'class_label'] = 'is_monomer_surface'
    # pose_df.to_csv(name_no_suffix + ".dat", index=False, sep=" ")


    return pose_df, atomic_depth if interface_mode else atomic_depth_monomer



scorefxn_sc = core.scoring.ScoreFunctionFactory.create_score_function("none")
scorefxn_sc.set_weight(core.scoring.hbond_bb_sc, 1)
scorefxn_sc.set_weight(core.scoring.hbond_sc, 1)

scorefxn_bb = core.scoring.ScoreFunctionFactory.create_score_function("none")
scorefxn_bb.set_weight(core.scoring.hbond_sr_bb, 1)
scorefxn_bb.set_weight(core.scoring.hbond_lr_bb, 1)


def find_his_dist_pairs(pose, rotsets, interface_mode=False, dist_cutoff_mode=-1):

    monomer_size = pose.conformation().chain_end(1)

    print("Locating polar atoms")
    rotamer_polar_atoms = [None]
    rotamer_polar_iatoms = [None]
    for irotamer in range(1, rotsets.nrotamers()+1):
        rotamer = rotsets.rotamer(irotamer)
        if rotamer.name1() in 'DERKH':
            these_atoms = []
            these_iatoms = []
            for iatom in range(rotamer.first_sidechain_atom(), rotamer.nheavyatoms()+1):
                if rotamer.atom_type(iatom).element() != 'C':
                    these_atoms.append(nup.from_vector(rotamer.xyz(iatom)))
                    these_iatoms.append(iatom)
            these_atoms = np.array(these_atoms)
            these_iatoms = np.array(these_iatoms)
        else:
            these_atoms = np.zeros((0,3), dtype=float)
            these_iatoms = np.zeros((0,), dtype=int)
        rotamer_polar_atoms.append(these_atoms)
        rotamer_polar_iatoms.append(these_iatoms)

    # store:
    # name1
    # irotamer
    # atom locs
    # all atom locs by name1
    atom_dict = defaultdict( lambda : defaultdict( lambda : dict(all_atoms=[], ind_rotamers=[])))
    for irotamer in range(1, rotsets.nrotamers()+1):
        rotamer = rotsets.rotamer(irotamer)
        seqpos = rotsets.res_for_rotamer(irotamer)
        these_atoms = rotamer_polar_atoms[irotamer]
        these_iatoms = rotamer_polar_iatoms[irotamer]
        name1 = rotamer.name1()
        atom_dict[seqpos][name1]['all_atoms'].append(these_atoms)
        atom_dict[seqpos][name1]['ind_rotamers'].append((irotamer, these_atoms, these_iatoms))

    for seqpos in atom_dict:
        for name1 in atom_dict[seqpos]:
            atom_dict[seqpos][name1]['atom_starts'] = np.r_[0, np.cumsum([len(x) for x in atom_dict[seqpos][name1]['all_atoms'][:-1]])]
            atom_dict[seqpos][name1]['all_atoms'] = np.concatenate(atom_dict[seqpos][name1]['all_atoms'])


    npose = nup.npose_from_pose(pose)
    Cas = nu.extract_atoms(npose, [nu.CA])

    Ca_dist = np.linalg.norm(Cas[:,None] - Cas[None,:], axis=-1)

    dist_cutoff = 25
    Ca_close = Ca_dist < dist_cutoff
    is_binder = np.zeros(pose.size(), dtype=bool)
    is_binder[:monomer_size] = True


    if interface_mode:
        chain_ok = is_binder[:,None] != is_binder[None,:]
    else:
        chain_ok = is_binder[:,None] & is_binder[None,:]


    ok_mask = Ca_close & chain_ok
    wh_ok = np.where(ok_mask)
    wh_ok_mask = wh_ok[0] < wh_ok[1]
    wh_ok = (wh_ok[0][wh_ok_mask] + 1, wh_ok[1][wh_ok_mask] + 1) # mask and to 1-indexed

    dist_cutoff2 = dist_cutoff_mode**2

    hb_pairs = {}
    print('Doing all by all rotamer distances')
    for seqpos, other_seqpos in zip(*wh_ok):

        for name1 in atom_dict[seqpos]:
            for other_name1 in atom_dict[other_seqpos]:

                if name1 != 'H' and other_name1 != 'H':
                    continue

                # No mutation lol
                if name1 == pose.residue(seqpos).name1() and other_name1 == pose.residue(other_seqpos).name1():
                    continue

                # this happens in interface mode
                if other_name1 not in ('H' + args.other_aas):
                    continue

                d2_map = np.sum(np.square(atom_dict[seqpos][name1]['all_atoms'][:,None] - atom_dict[other_seqpos][other_name1]['all_atoms'][None,:]), axis=-1)

                if d2_map.min() > dist_cutoff2:
                    continue

                _iatom, other__iatom = np.unravel_index(np.argmin(d2_map), d2_map.shape)

                local_irot = np.searchsorted(atom_dict[seqpos][name1]['atom_starts'], _iatom, 'right')-1
                other_local_irot = np.searchsorted(atom_dict[other_seqpos][other_name1]['atom_starts'], other__iatom, 'right')-1

                irotamer, our_atoms, our_iatoms = atom_dict[seqpos][name1]['ind_rotamers'][local_irot]
                other_irotamer, other_our_atoms, other_our_iatoms = atom_dict[other_seqpos][other_name1]['ind_rotamers'][other_local_irot]

                ind_d2 = np.sum(np.square(our_atoms[:,None] - other_our_atoms[None,:]), axis=-1)
                closest = ind_d2.min()
                assert np.isclose(d2_map.min(), ind_d2.min())

                idx_atom, other_idx_atom = np.unravel_index(np.argmin(ind_d2), ind_d2.shape)
                iatom = our_iatoms[idx_atom]
                other_iatom = other_our_iatoms[other_idx_atom]

                rotamer = rotsets.rotamer(irotamer)
                other_rotamer = rotsets.rotamer(other_irotamer)
                assert np.isclose(rotamer.xyz(iatom).distance(other_rotamer.xyz(other_iatom)), np.sqrt(closest))

                # storage

                name1seqpos = name1 + str(seqpos)
                other_name1seqpos = other_name1 + str(other_seqpos)

                if name1seqpos < other_name1seqpos:
                    name1seqpos1, name1seqpos2 = name1seqpos, other_name1seqpos
                    rotamer1, rotamer2 = rotamer, other_rotamer
                    atom1, atom2 = iatom, other_iatom
                else:
                    name1seqpos2, name1seqpos1 = name1seqpos, other_name1seqpos
                    rotamer2, rotamer1 = rotamer, other_rotamer
                    atom2, atom1 = iatom, other_iatom

                key = (name1seqpos1, name1seqpos2)

                store = True
                e = np.sqrt(closest)
                if key in hb_pairs:
                    last_e = hb_pairs[key][0]
                    if last_e < e:
                        store = False

                if store:
                    hb_pairs[key] = [e, rotamer1, rotamer2, iatom, other_iatom]


    print("Found %i pairs"%(len(hb_pairs)))

    return hb_pairs


def find_his_hbond_pairs(pose, interface_mode=False, two_sided_design=False, dist_cutoff_mode=-1):

    monomer_size = pose.conformation().chain_end(1)

    tf = core.pack.task.TaskFactory()

    if not two_sided_design:
        repack_aa = core.pack.task.operation.RestrictToRepackingRLT()
        subset = chainB.apply(pose)
        tf.push_back( core.pack.task.operation.OperateOnResidueSubset( repack_aa, subset ) )
    else:
        restrict_aa = core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
        restrict_aa.aas_to_keep( 'H' + args.other_aas )
        subset = chainB.apply(pose)
        # Disulfide cysteines can't be designed and would end up with no rotamers (which breaks the hbond graph),
        #  so only repack them
        disulfide_subset = utility.vector1_bool(pose.size())
        for seqpos in range(1, pose.size()+1):
            if pose.residue(seqpos).type().is_disulfide_bonded():
                disulfide_subset[seqpos] = True
                subset[seqpos] = False
        tf.push_back( core.pack.task.operation.OperateOnResidueSubset( core.pack.task.operation.RestrictToRepackingRLT(), disulfide_subset ) )
        tf.push_back( core.pack.task.operation.OperateOnResidueSubset( restrict_aa, subset ) )


    restrict_aa = core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
    restrict_aa.aas_to_keep( 'H' + args.other_aas )
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
    print("Building rotamers")
    rotsets.build_rotamers( pose, scorefxn, graph )

    if dist_cutoff_mode > 0:
        return find_his_dist_pairs(pose, rotsets, interface_mode=interface_mode, dist_cutoff_mode=dist_cutoff_mode)


    complete_rotamer_sets = core.pack.rotamer_set.RotamerSets()
    position_had_rotset = utility.vector1_bool()

    print('Building hbgraph')
    hb_graph = core.pack.hbonds.hbond_graph_from_partial_rotsets(pose, rotsets, scorefxn_sc, scorefxn_bb, complete_rotamer_sets, position_had_rotset, -0.01)

    print("Processing")

    assert np.all(position_had_rotset)
    assert rotsets.nrotamers() == hb_graph.num_nodes()

    hb_pairs = {}

    for ihbnode in range(1, hb_graph.num_nodes()+1):

        hbnode = hb_graph.get_node(ihbnode)
        seqpos = rotsets.res_for_rotamer(ihbnode)
        rotamer = rotsets.rotamer(ihbnode)
        name1 = rotamer.name1()


        if seqpos > monomer_size:
            continue

        name1seqpos = name1 + str(seqpos)

        it = hbnode.edge_list_begin( hb_graph )
        while it.valid():
            edge = _hbedge_from_lowmem(hb_graph, it.dereference(), ihbnode)
            it.pre_increment()

            we_are_first = edge.get_first_node_ind() == ihbnode

            other_node = edge.get_second_node_ind() if we_are_first else edge.get_first_node_ind()
            assert other_node != ihbnode
            other_seqpos = rotsets.res_for_rotamer(other_node)
            other_rotamer = rotsets.rotamer(other_node)
            other_name1 = other_rotamer.name1()
            other_name1seqpos = other_name1 + str(other_seqpos)

            # At least one of them has to be H
            if name1 != 'H' and other_name1 != 'H':
                continue

            # No mutation lol
            if name1 == pose.residue(seqpos).name1() and other_name1 == pose.residue(other_seqpos).name1():
                continue

            # this happens in interface mode
            if other_name1 not in ('H' + args.other_aas):
                continue

            if interface_mode:
                if other_seqpos <= monomer_size:
                    continue
            else:
                if other_seqpos > monomer_size:
                    continue

            for hbond in edge.hbonds():

                we_are_donor = not (we_are_first ^ hbond.first_node_is_donor())

                if we_are_donor:
                    continue

                our_iatom = hbond.local_atom_id_A()
                their_iatom = hbond.local_atom_id_D()

                if rotamer.atom_is_backbone(our_iatom):
                    continue

                if other_rotamer.atom_is_backbone(their_iatom):
                    continue

                if name1seqpos < other_name1seqpos:
                    name1seqpos1, name1seqpos2 = name1seqpos, other_name1seqpos
                    rotamer1, rotamer2 = rotamer, other_rotamer
                    atom1, atom2 = our_iatom, their_iatom
                else:
                    name1seqpos2, name1seqpos1 = name1seqpos, other_name1seqpos
                    rotamer2, rotamer1 = rotamer, other_rotamer
                    atom2, atom1 = our_iatom, their_iatom

                key = (name1seqpos1, name1seqpos2)

                store = True
                e = hbond.score()
                if key in hb_pairs:
                    last_e = hb_pairs[key][0]
                    if last_e < e:
                        store = False

                if store:
                    hb_pairs[key] = [e, rotamer1, rotamer2, atom1, atom2]

    print("Found %i pairs"%(len(hb_pairs)))

    return hb_pairs






def worst_possible_asp(pose, name_no_suffix, out_score_map, out_string_map, suffix):

    class_df, atomic_depth_monomer = classify_binder_positions(pose, interface_mode=args.interface_mode)
    type_set = pose.residue(1).type().atom_type_set()

    monomer = pose.split_by_chain()[1]
    dssp = better_dssp3(monomer)
    ss_elems = get_ss_elements2(dssp)
    monomer_size = monomer.size()

    if args.interface_mode:
        target= pose.split_by_chain()[2]
        target_dssp = better_dssp3(target)
        target_ss_elems = [(tp, x+monomer_size, y+monomer_size) for tp, x, y in get_ss_elements2(target_dssp)]
        dssp += target_dssp[1:]
        ss_elems += target_ss_elems


    on_element = np.full((len(dssp),), -1000, dtype=int)
    for ielem, (letter, start, end) in enumerate(ss_elems):
        on_element[start:end+1] = ielem

    elem_is_interface = np.zeros(len(ss_elems), dtype=bool)
    for ielem, (letter, start, end) in enumerate(ss_elems):
        rows = class_df[(class_df['seqpos_int'] >= start) & (class_df['seqpos_int'] <= end)]
        any_interface = (rows['is_interface_core'] | rows['is_interface_boundary']).any()
        elem_is_interface[ielem] = any_interface


    his_pairs = find_his_hbond_pairs(pose, interface_mode=args.interface_mode, two_sided_design=args.two_sided_design, 
                                                                                        dist_cutoff_mode=args.dist_cutoff_mode)

    pose_copy = pose.clone()

    out_stuff = []

    only_allow_set = None
    if args.only_allow_positions:
        only_allow_set = set()
        for value in args.only_allow_positions.split(','):
            value = int(value)
            assert value > 0, '--only_allow_positions must be >= 0'
            assert value <= monomer_size, f'--only_allow_positions are 1-indexed. Your monomer has max size {monomer_size} but you input {value}'
            only_allow_set.add(value)

    for (name1seqpos1, name1seqpos2), (e, rotamer1, rotamer2, iatom1, iatom2) in his_pairs.items():

        seqpos1 = int(name1seqpos1[1:])
        seqpos2 = int(name1seqpos2[1:])
        name1_1 = name1seqpos1[:1]
        name1_2 = name1seqpos2[:1]

        if only_allow_set:
            if seqpos1 not in only_allow_set:
                continue
            if seqpos2 not in only_allow_set:
                continue

        out_score_map = std.map_std_string_double()
        out_string_map = std.map_std_string_std_string()
        tag = name_no_suffix + f'_{seqpos1}{name1_1}_{seqpos2}{name1_2}'

        pose = pose_copy.clone()

        row1 = class_df[class_df['seqpos'] == str(seqpos1)].iloc[0]
        row2 = class_df[class_df['seqpos'] == str(seqpos2)].iloc[0]

        class1 = row1['class_label'].replace('is_', '')
        class2 = row2['class_label'].replace('is_', '')

        energy_label = 'min_distance' if args.dist_cutoff_mode > 0 else 'energy'
        out_score_map[energy_label] = e

        max_sasa = max(row1['sc_sasa'], row2['sc_sasa'])
        out_score_map['max_sc_sasa'] = max_sasa

        class_a, class_b = list(sorted([class1, class2]))

        out_string_map['class'] = class_a + '-' + class_b

        pose.replace_residue(seqpos1, rotamer1.clone(), True)
        pose.replace_residue(seqpos2, rotamer2.clone(), True)
        pose.pdb_info().add_reslabel(seqpos1, 'internal_HIS')
        pose.pdb_info().add_reslabel(seqpos2, 'internal_HIS')

        atom1 = rotamer1.atom(iatom1)
        atom2 = rotamer2.atom(iatom2)

        depth1 = atomic_depth_monomer.calcdepth(atom1, type_set)
        depth2 = atomic_depth_monomer.calcdepth(atom2, type_set)

        out_score_map['hbond_depth'] = depth1 + depth2

        dssp1 = dssp[seqpos1]
        dssp2 = dssp[seqpos2]

        out_string_map['dssp'] = ''.join(list(sorted([dssp1, dssp2])))

        elem1 = on_element[seqpos1]
        elem2 = on_element[seqpos2]

        out_score_map['separate_ss_elems'] = elem1 != elem2
        out_score_map['n_interface_elems'] = int(elem_is_interface[elem1]) + int(elem_is_interface[elem2])

        out_score_map['seqpos1'] = seqpos1
        out_score_map['seqpos2'] = seqpos2
        out_string_map['aa1'] = name1_1
        out_string_map['aa2'] = name1_2


        print(f"{rotamer1.name3()}-{rotamer2.name3()} {seqpos1}-{seqpos2}, dssp: {out_string_map['dssp']} n_interface: {out_score_map['n_interface_elems']}"
                f" separate: {out_score_map['separate_ss_elems']} class: {out_string_map['class']} depth: {out_score_map['hbond_depth']:.1f}"
                f" max_sasa: {out_score_map['max_sc_sasa']:.1f} {energy_label}: {out_score_map[energy_label]:.1f}")

        out_stuff.append([pose, tag, out_score_map, out_string_map])

    return out_stuff



############### BEGIN MAIN FUNCTION ###########################

if ( silent != "" ):
    sfd_in = rosetta.core.io.silent.SilentFileData(rosetta.core.io.silent.SilentFileOptions())
    sfd_in.read_file(silent)

    pdbs = list(sfd_in.tags())

    sfd_out = core.io.silent.SilentFileData( "out.silent", False, False, "binary", core.io.silent.SilentFileOptions())



# ckpt = 0
# if (os.path.exists("ckpt")):
#     with open("ckpt") as f:
#         try:
#             ckpt = int(f.read().strip())
#         except:
#             pass


num = -1
for ipdb, pdb in enumerate(pdbs):
    t0 = time.time()
    print("Attempting pose: " + pdb)


    # if ( ckpt > ipdb ):
    #     print("Checkpoint: Already finished: " + pdb)
    #     continue
    # with open("ckpt", "w") as f:
    #     f.write(str(ipdb))

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


        out_pose = worst_possible_asp(pose, name_no_suffix, score_map, string_map, "")

        to_iterate = [(out_pose, name_no_suffix, score_map, string_map)]

        if ( isinstance(out_pose, list) ):
            to_iterate = out_pose

        for pose, this_name_no_suffix, score_map, string_map in to_iterate:

            core.io.raw_data.ScoreMap.add_arbitrary_score_data_from_pose( pose, score_map)
            core.io.raw_data.ScoreMap.add_arbitrary_string_data_from_pose( pose, string_map)
            sfd.write_pose( pose, score_map, this_name_no_suffix, string_map)

            if ( silent == "" ):
                pose.dump_pdb(this_name_no_suffix + ".pdb")
            else:
                silent_name = "out.silent"
                sfd_out = core.io.silent.SilentFileData( silent_name, False, False, "binary", core.io.silent.SilentFileOptions())
                struct = sfd_out.create_SilentStructOP()
                struct.fill_struct(pose, this_name_no_suffix)
                for score in score_map:
                    struct.add_energy(score, score_map[score])
                for string in string_map:
                    struct.add_string_value(string, string_map[string])
                sfd_out.add_structure(struct)
                sfd_out.write_all(silent_name, False)


        seconds = int(time.time() - t0)

        print("protocols.jd2.JobDistributor: " + name_no_suffix + " reported success in %i seconds"%seconds)

    # except Exception as e:
    #     print("Error!!!")
    #     print(e)



# if ( silent != "" ):
#     sfd_out.write_all("out.silent", False)




















