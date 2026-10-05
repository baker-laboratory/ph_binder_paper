#!/usr/bin/env python
from __future__ import division

# This program accepts arguments like this:

#./remove_superfluous_trp.py pdb1.pdb pdb2.pdb pdb3.pdb
# or
#./remove_superfluous_trp.py -in:file:silent my.silent

import os
import sys
import math

import distutils.spawn
import os
import sys
#sys.path.append(os.path.dirname(distutils.spawn.find_executable("silent_tools.py")))
#import silent_tools

from pyrosetta import *
from pyrosetta.rosetta import *

sys.path.append("/home/bcov/sc/random/npose")
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

# import pyRMSD.RMSDCalculator

init("-beta_nov16 -in:file:silent_struct_type binary -keep_input_scores false -mute all"
    " -holes:dalphaball /work/tlinsky/Rosetta/main/source/external/DAlpahBall/DAlphaBall.macgcc"
    )




parser = argparse.ArgumentParser()
parser.add_argument("-in:file:silent", type=str, default="")
parser.add_argument("pdbs", type=str, nargs="*")
parser.add_argument("--list_of_mpnn_prob_npz", type=str, default='')
parser.add_argument("--modes", type=str, default='regular,force_touching')
parser.add_argument("--two_sided_design", action='store_true')

args = parser.parse_args(sys.argv[1:])

mpnn_npz_paths = {}
with open(args.list_of_mpnn_prob_npz) as f:
    for line in f:
        line = line.strip()
        if len(line) == 0:
            continue
        mpnn_npz_paths[os.path.basename(line).replace('.npz', '')] = line



pdbs = args.pdbs
silent = args.__getattribute__("in:file:silent")



def fix_scorefxn(sfxn, allow_double_bb=False):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    sfxn.set_energy_method_options(opts)

scorefxn = get_fa_scorefxn()
sfxn_design = scorefxn.clone()
sfxn_design.set_weight(core.scoring.res_type_constraint, 1)
scorefxn_fa_atr = core.scoring.ScoreFunctionFactory.create_score_function("none")
scorefxn_fa_atr.set_weight(core.scoring.fa_atr, 1)
scorefxn_vdw = core.scoring.ScoreFunctionFactory.create_score_function("none")
scorefxn_vdw.set_weight(core.scoring.fa_atr, 1)
scorefxn_vdw.set_weight(core.scoring.fa_rep, 0.55)

sfxn_cart = core.scoring.ScoreFunctionFactory.create_score_function("beta_nov16_cart")

fix_scorefxn(scorefxn)

all_weights = scorefxn.get_nonzero_weighted_scoretypes()
lr_terms = core.scoring.EMapVector()
lr_terms.assign(scorefxn.weights())
for weight in all_weights:
    if ( "dslf"  not in core.scoring.name_from_score_type(weight) and "rama_prepro" not in core.scoring.name_from_score_type(weight)):
        lr_terms.set(weight, 0)

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



# Newer versions of pyrosetta can't load the FSP like that
#xml = '''
#<MOVERS>
#    <FavorSequenceProfile name="upweight_native" scaling="global" weight="1" chain="1" use_starting="true" matrix="MATCH" />
#</MOVERS>

#'''

#objs = protocols.rosetta_scripts.XmlObjects.create_from_string(xml)
# upweight_native = objs.get_mover('upweight_native')


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



def grab_energies(pose, e_onebody, e_twobody, e_twobody_gly):

    mut_pos = 0

    if ( mut_pos is None ):
        print("Bad seq")
        return None

    n = 0

    whole_pose = core.select.residue_selector.TrueResidueSelector().apply(pose)
    gly_pose = pose.clone()
    protocols.toolbox.pose_manipulation.repack_these_residues(whole_pose, gly_pose, scorefxn, False, "G")

    scorefxn(pose)
    scorefxn(gly_pose)

    for seqpos in range(1, pose.size()+1):
        if ( seqpos != mut_pos ):
            # e_onebody[n, seqpos] = pose.energies().residue_total_energies(seqpos).dot(oneb_weights)
            e_onebody[n, seqpos] = pose.energies().onebody_energies(seqpos).dot(scorefxn.weights())
            e_onebody[n, seqpos] += pose.energies().residue_total_energies(seqpos).dot(lr_terms)
            # e_onebody[n, seqpos] = pose.energies().residue_total_energies(seqpos).dot(scorefxn.weights())

    for seqpos1 in range(1, pose.size()+1):
        if ( seqpos1 == mut_pos ):
            continue
        for seqpos2 in range(1, pose.size()+1):
            if ( seqpos2 <= seqpos1 ):
                continue
            edge = pose.energies().energy_graph().find_edge(seqpos1, seqpos2)
            if ( edge ):
                emap = edge.fill_energy_map()
                e_twobody[n, seqpos1, seqpos2] = emap.dot(scorefxn.weights())
            else:
                e_twobody[n, seqpos1, seqpos2] = 0


            edge = gly_pose.energies().energy_graph().find_edge(seqpos1, seqpos2)
            if ( edge ):
                emap = edge.fill_energy_map()
                e_twobody_gly[n, seqpos1, seqpos2] = emap.dot(scorefxn.weights())
            else:
                e_twobody_gly[n, seqpos1, seqpos2] = 0

    if ( mut_pos == 0 ):
        assert(np.isclose(np.nan_to_num(e_twobody[n]).sum() + np.nan_to_num(e_onebody[n]).sum(), scorefxn(pose)))


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



def do_a_ppi_fast_design(pose, prefix, suffix, to_design, min_mask, monomer_size, native_weight, mpnn_probs, alphabet, mpnn_cutoff=-2.5, fast=False, two_sided_design=False):

    pose = pose.clone()


    mpnn_ok = mpnn_probs > mpnn_cutoff

    score_map = std.map_std_string_double()
    string_map = std.map_std_string_std_string()


    locked = ~to_design
    locked_subset = utility.vector1_bool( pose.size() )
    for i in range(1, pose.size()+1):
        if i <= monomer_size or two_sided_design:
            locked_subset[i] = locked[i-1]
        else:
            locked_subset[i] = True

    tf = core.pack.task.TaskFactory()
    tf.push_back( core.pack.task.operation.OperateOnResidueSubset( core.pack.task.operation.PreventRepackingRLT(),
                  locked_subset ) )

    for seqpos0 in np.where(to_design)[0]:
        seqpos = seqpos0 + 1

        ok_mask = mpnn_ok[seqpos0]
        ok_letters = set([alphabet[x] for x in np.where(ok_mask)[0]])
        ok_letters.add(pose.sequence()[seqpos0])

        if len(ok_letters) < 3:
            arg_weights = np.argsort(-mpnn_probs[seqpos0])
            for iletter in arg_weights[:3]:
                ok_letters.add(alphabet[iletter])

        ok_letters = ''.join(list(ok_letters))

        restrict_aa = core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
        restrict_aa.aas_to_keep( ok_letters )

        subset = utility.vector1_bool( pose.size() )
        subset[seqpos] = True
        tf.push_back( core.pack.task.operation.OperateOnResidueSubset( restrict_aa, subset ) )


    tf.push_back( core.pack.task.operation.IncludeCurrent() )

    ex1_ex2 = core.pack.task.operation.ExtraRotamersGeneric()
    ex1_ex2.ex1( True )
    ex1_ex2.ex2aro( True )
    tf.push_back(ex1_ex2)

    packer = protocols.minimization_packing.PackRotamersMover()
    packer.score_function( sfxn_design )
    packer.task_factory(tf)

    # Build the mover directly, bypassing XML parse
    upweight_native = protocols.simple_moves.FavorSequenceProfile()
    upweight_native.set_weight(1.0)
    upweight_native.set_scaling("global")

    # Capture sequence from starting pose (equivalent to use_starting=true)
    seq = core.sequence.Sequence(pose)
    upweight_native.set_sequence(seq, "MATCH")    

    upweight_native.set_weight(native_weight)
    upweight_native.apply(pose)


    mm = core.kinematics.MoveMap()
    for seqpos in range(1, pose.size()+1):
        seqpos0 = seqpos-1
        if (seqpos0 < monomer_size or two_sided_design) and min_mask[seqpos0]:
            mm.set_bb(seqpos, 1)
            mm.set_chi(seqpos, 1)
        else:
            mm.set_bb(seqpos, 0)
            mm.set_chi(seqpos, 0)

    min_mover = protocols.minimization_packing.MinMover()
    min_mover.set_movemap( mm )
    min_mover.tolerance( 0.01 )
    min_mover.score_function( sfxn_cart )
    min_mover.min_type( "lbfgs_armijo_nonmonotone" )
    min_mover.cartesian(True)
    min_mover.min_options().max_iter(200)


    jack_maguire_weights = [0.079, 0.295, 0.577, 1]
    if fast:
        jack_maguire_weights = [0.35, 1] # make the 0.35 so we don't trigger the dumps

    out_stuff = []
    for weight in jack_maguire_weights:
        sfxn_design.set_weight(core.scoring.fa_rep, 0.55*weight)
        packer.score_function(sfxn_design)

        sfxn_cart.set_weight(core.scoring.fa_rep, 0.55*weight)
        min_mover.score_function(sfxn_cart)

        if np.abs(weight - 0.295) < 0.01:
            out_stuff.append([pose.clone(), prefix + '_early3' + suffix, score_map, string_map])

        if np.abs(weight - 0.577) < 0.01:
            out_stuff.append([pose.clone(), prefix + '_early2' + suffix, score_map, string_map])

        if weight == 1:
            out_stuff.append([pose.clone(), prefix + '_early1' + suffix, score_map, string_map])

        print("Pack:", weight)
        packer.apply(pose)
        if weight < 1:
            print('Min:', weight)
            min_mover.apply(pose)


    out_stuff.append([pose, prefix + suffix, score_map, string_map])

    return out_stuff



def which_aa_can_touch_his(pose, his_positions, ok_letters, touch_cutoff=-0.5, max_allowed=20):

    pose = pose.clone()

    design_subset = utility.vector1_bool( pose.size() )
    for seqpos in ok_letters:
        design_subset[seqpos] = True

    # Mutate the design residues to GLY so we can look at them one by one
    protocols.toolbox.pose_manipulation.repack_these_residues(design_subset, pose, scorefxn, False, 'G')


    best_scores_at_seqpos = {}
    for test_seqpos, letters in ok_letters.items():

        had_g = 'G' in letters
        tf = core.pack.task.TaskFactory()

        locked_subset = utility.vector1_bool( pose.size() )
        for seqpos in range(1, pose.size()+1):
            locked_subset[seqpos] = True
        for his_seqpos in his_positions:

            restrict_aa = core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
            restrict_aa.aas_to_keep( pose.residue(his_seqpos).name1() )

            subset = utility.vector1_bool( pose.size() )
            subset[his_seqpos] = True
            tf.push_back( core.pack.task.operation.OperateOnResidueSubset( restrict_aa, subset ) )

            locked_subset[his_seqpos] = False


        restrict_aa = core.pack.task.operation.RestrictAbsentCanonicalAASRLT()
        restrict_aa.aas_to_keep( letters + 'G')

        subset = utility.vector1_bool( pose.size() )
        subset[test_seqpos] = True
        tf.push_back( core.pack.task.operation.OperateOnResidueSubset( restrict_aa, subset ) )

        locked_subset[test_seqpos] = False

        tf.push_back( core.pack.task.operation.OperateOnResidueSubset( core.pack.task.operation.PreventRepackingRLT(),
                      locked_subset ) )

        tf.push_back( core.pack.task.operation.IncludeCurrent() )



        rotamer_sets = core.pack.rotamer_set.RotamerSets()

        ig = core.pack.pack_rotamers_setup(pose, scorefxn_vdw, tf.create_task_and_apply_taskoperations( pose ), rotamer_sets)


        his_og_rots = []
        for his_seqpos in his_positions:
            res = pose.residue(his_seqpos)
            assert rotamer_sets.has_rotamer_set_for_residue(his_seqpos)
            rotset = rotamer_sets.rotamer_set_for_residue(his_seqpos)
            found_irot = None
            for irot in range(1, rotset.num_rotamers()+1):
                matches = True
                rotamer = rotset.rotamer(irot)
                assert rotamer.name1() == res.name1(), f'{rotamer.name1()} == {res.name1()}'
                for ichi in range(1, res.nchi()+1):
                    if np.abs( rotamer.chi(ichi) - res.chi(ichi) )%360 > 1:
                        matches = False
                        break
                if matches:
                    found_irot = irot
                    break
            assert found_irot is not None, "Rotamer not found?"

            his_og_rots.append(found_irot)

        molten_2_seqpos = rotamer_sets.moltenres_2_resid_vector()

        his_moltenres = [list(molten_2_seqpos).index(x)+1 for x in his_positions]

        test_moltenres = list(molten_2_seqpos).index(test_seqpos)+1


        test_rotset = rotamer_sets.rotamer_set_for_residue(test_seqpos)

        letter_best_score = {letter:0 for letter in letters}

        for ihis in range(len(his_positions)):
            his_molt = his_moltenres[ihis]
            his_irot = his_og_rots[ihis]

            edge = ig.find_edge(test_moltenres, his_molt)
            if edge is None:
                continue

            letter_best_score_inner = {letter:0 for letter in letters}
            letter_best_score_inner['G'] = 0

            test_is_first = edge.get_first_node_ind() == test_moltenres
            if test_is_first:
                assert edge.get_first_node_ind() == test_moltenres
                assert edge.get_second_node_ind() == his_molt
            else:
                assert edge.get_second_node_ind() == test_moltenres
                assert edge.get_first_node_ind() == his_molt

            for irot in range(1, test_rotset.num_rotamers()+1):
                test_rotamer = test_rotset.rotamer(irot)
                letter = test_rotamer.name1()

                if test_is_first:
                    energy = edge.get_two_body_energy(irot, his_irot)
                else:
                    energy = edge.get_two_body_energy(his_irot, irot)

                if energy < touch_cutoff or letter == 'G':
                    letter_best_score_inner[letter] = min(letter_best_score_inner[letter], energy)

            G_energy = letter_best_score_inner['G']
            for letter, score in letter_best_score_inner.items():
                better_than_G = score - G_energy
                if better_than_G < touch_cutoff:
                    letter_best_score[letter] = better_than_G


        best_scores_at_seqpos[test_seqpos] = letter_best_score



    for cutoff in np.linspace(touch_cutoff, touch_cutoff-100, 1000):
        N_letters = 0
        for seqpos in best_scores_at_seqpos:
            for letter, score in best_scores_at_seqpos[seqpos].items():
                if score < cutoff:
                    N_letters += 1

        if N_letters <= max_allowed:
            break


    interacting_letters = {}
    for seqpos in best_scores_at_seqpos:
        new_ok_letters = ''
        for letter, score in best_scores_at_seqpos[seqpos].items():
            if score < cutoff:
                new_ok_letters += letter
        if len(new_ok_letters) > 0:
            interacting_letters[seqpos] = new_ok_letters





    print('Found these residues to interact with HIS')
    for seqpos, letters in interacting_letters.items():
        print(seqpos, letters)

        # import IPython
        # IPython.embed()
    return interacting_letters



the_locals = None

def worst_possible_asp(pose, name_no_suffix, out_score_map, out_string_map, suffix):

    monomer_size = pose.conformation().chain_end(1)
    look_size = pose.size() if args.two_sided_design else monomer_size

    class_df, _ = classify_binder_positions(pose, interface_mode=args.two_sided_design)

    npose = nup.npose_from_pose(pose)

    his_positions = []
    backup_positions = []
    for seqpos in range(1, look_size+1):
        reslabels = pose.pdb_info().get_reslabels(seqpos)
        if 'internal_HIS' in list(reslabels):
            his_positions.append(seqpos)
        for label in reslabels:
            if 'backup_' in label:
                backup_positions.append(seqpos)

    assert len(his_positions) == 2, his_positions

    his_positions = np.array(his_positions, dtype=int)
    backup_positions = np.array(backup_positions, dtype=int)
    network_positions = np.concatenate((his_positions, backup_positions))

    seqsize = pose.size()

    e_onebody = np.zeros((1+seqsize))
    e_onebody[:] = np.nan


    e_twobody = np.zeros((1+seqsize, 1+seqsize))
    e_twobody[:] = np.nan

    e_twobody_gly = np.zeros((1+seqsize, 1+seqsize))
    e_twobody_gly[:] = np.nan


    # Only the e_array's are 1-indexed. We switch everything to 0 index after
    grab_energies(pose, e_onebody[None], e_twobody[None], e_twobody_gly[None])


    e_twobody -= e_twobody_gly

    per_res_ddg = np.concatenate((e_twobody[1:monomer_size+1,monomer_size+1:].sum(axis=-1), (e_twobody[1:monomer_size+1,monomer_size+1:].sum(axis=0))))

    is_interface = (np.abs(per_res_ddg) > 0.25)[:look_size]

    clashing_with_his = (e_twobody[1:look_size+1][:,network_positions] > 0.5).any(axis=-1)


    dist_2_his = np.zeros(look_size)
    his_atoms = []
    for seqpos in network_positions:
        res = pose.residue(seqpos)
        for iatom in range(res.first_sidechain_atom(), res.nheavyatoms()+1):
            his_atoms.append(nup.from_vector(res.xyz(iatom)))
    his_atoms = np.array(his_atoms)

    for seqpos in range(1, look_size+1):
        res = pose.residue(seqpos)
        atoms = [nup.from_vector(res.nbr_atom_xyz())]
        for iatom in range(res.first_sidechain_atom(), res.nheavyatoms()+1):
            atoms.append(nup.from_vector(res.xyz(iatom)))

        atoms = np.array(atoms)
        dist_2_his[seqpos-1] = np.linalg.norm(atoms[:,None] - his_atoms[None,:], axis=-1).min()


    Cb = nu.extract_atoms(npose, [nu.CB])
    Cb_dist_his = np.linalg.norm( Cb[:,None] - Cb[None,network_positions-1], axis=-1).min(axis=-1)


    mpnn_file = mpnn_npz_paths[name_no_suffix]
    asdf = np.load(mpnn_file)
    mpnn_probs = asdf['odds'][0]
    alphabet = str(asdf['alphabet'])

    mpnn_probs[:,alphabet.index('C')] = -1000
    mpnn_probs[:,alphabet.index('X')] = -1000


    mpnn_cutoff = -2.5
    mpnn_ok = mpnn_probs > mpnn_cutoff

    seq = pose.sequence()
    is_ala = np.array([(x in 'AG') for x in seq])[:look_size]

    pose_copy = pose.clone()

    out_stuff = []

    if 'regular' in args.modes.split(','):
        for his_dist in [5, 5.5, 6, 6.5]:
            for native_weight in [-1, -2, -3]:

                pose = pose_copy.clone()

                to_design = (clashing_with_his | (dist_2_his < his_dist) | ((dist_2_his < his_dist + 2.5) & is_ala) ) & ~is_interface
                to_design[network_positions-1] = False

                min_mask = to_design | (Cb_dist_his[:look_size] < 9)
                # And then 2 on each side of all designed residues
                min_mask[1:] |= to_design[:-1]
                min_mask[2:] |= to_design[:-2]
                min_mask[:-1] |= to_design[1:]
                min_mask[:-2] |= to_design[2:]

                print("HIS:", '+'.join([str(x) for x in his_positions]))
                print("Network:", '+'.join([str(x) for x in network_positions]))
                print("Clash:", '+'.join([str(x) for x in np.where(clashing_with_his)[0]+1]))
                print("Design:", '+'.join([str(x) for x in np.where(to_design)[0]+1]))
                print("Minimize:", '+'.join([str(x) for x in np.where(min_mask)[0]+1]))

                prefix = name_no_suffix + f'_dist{his_dist}_nat{-native_weight}'
                suffix = '_v4'

                out_stuff += do_a_ppi_fast_design(pose, prefix, suffix, to_design, min_mask, monomer_size, native_weight, mpnn_probs, alphabet, mpnn_cutoff=-2.5, two_sided_design=args.two_sided_design)
    
    if 'force_touching' in args.modes.split(','):

        skip_positions = []
        for seqpos in range(1, look_size+1):
            if class_df[class_df['seqpos'] == str(seqpos)].iloc[0]['is_monomer_surface']:
                skip_positions.append(seqpos)

        his_dist = 6
        native_weight = -2

        to_design = (clashing_with_his | (dist_2_his < his_dist) | ((dist_2_his < his_dist + 2.5) & is_ala) ) & ~is_interface
        to_design[network_positions-1] = False

        min_mask = to_design | (Cb_dist_his[:look_size] < 9)
        # And then 2 on each side of all designed residues
        min_mask[1:] |= to_design[:-1]
        min_mask[2:] |= to_design[:-2]
        min_mask[:-1] |= to_design[1:]
        min_mask[:-2] |= to_design[2:]

        ok_letters_dict = {}
        for seqpos0 in np.where(to_design)[0]:
            seqpos = seqpos0 + 1
            if seqpos in skip_positions:
                continue

            ok_mask = mpnn_ok[seqpos0]
            ok_letters = set([alphabet[x] for x in np.where(ok_mask)[0]])
            ok_letters.add(pose.sequence()[seqpos0])

            if len(ok_letters) < 3:
                arg_weights = np.argsort(-mpnn_probs[seqpos0])
                for iletter in arg_weights[:3]:
                    ok_letters.add(alphabet[iletter])

            ok_letters = ''.join(list(ok_letters))
            ok_letters_dict[seqpos] = ok_letters

        new_letters_dict = which_aa_can_touch_his(pose, his_positions, ok_letters_dict)

        trials = []
        for seqpos, letters in new_letters_dict.items():
            for letter in letters:
                trials.append((seqpos, letter))


        for trial in trials:
            test_seqpos, test_letter = trial

            prefix = name_no_suffix + f'_force{test_seqpos}{test_letter}'
            suffix = '_v4'
            pose = pose_copy.clone()

            to_design = to_design.copy()
            min_mask = min_mask.copy()

            to_design[test_seqpos-1] = False

            print("Force", f'{test_seqpos}{test_letter}')

            protocols.toolbox.pose_manipulation.repack_this_residue(test_seqpos, pose, scorefxn, False, test_letter)

            out_stuff += do_a_ppi_fast_design(pose, prefix, suffix, to_design, min_mask, monomer_size, native_weight, mpnn_probs, alphabet, mpnn_cutoff=-2.5, fast=True, two_sided_design=args.two_sided_design)



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




















