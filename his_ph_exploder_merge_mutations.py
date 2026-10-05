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

import pyrosetta
import pyrosetta.rosetta

def pyro():
    return pyrosetta
def ros():
    return pyrosetta.rosetta

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

init("-beta_nov16 -in:file:silent_struct_type binary -keep_input_scores false -mute all -ex1 -ex2"
    " -holes:dalphaball /work/tlinsky/Rosetta/main/source/external/DAlpahBall/DAlphaBall.macgcc"
    )




parser = argparse.ArgumentParser()
parser.add_argument("-in:file:silent", type=str, default="")
parser.add_argument("pdbs", type=str, nargs="*")
parser.add_argument("--original_pdb", type=str, default='')


args = parser.parse_args(sys.argv[1:])

pdbs = args.pdbs
silent = args.__getattribute__("in:file:silent")




og_pose = pose_from_file(args.original_pdb)


scorefxn = get_fa_scorefxn()
scorefxn_fa_atr = core.scoring.ScoreFunctionFactory.create_score_function("none")
scorefxn_fa_atr.set_weight(core.scoring.fa_atr, 1)

scorefxn_none = core.scoring.ScoreFunctionFactory.create_score_function("none")

def fix_scorefxn(sfxn, allow_double_bb=False):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    sfxn.set_energy_method_options(opts)


scorefxn_hbset = scorefxn.clone()
fix_scorefxn(scorefxn_hbset, True)


fix_scorefxn(scorefxn)
all_weights = scorefxn.get_nonzero_weighted_scoretypes()
lr_terms = core.scoring.EMapVector()
lr_terms.assign(scorefxn.weights())
for weight in all_weights:
    if ( "dslf"  not in core.scoring.name_from_score_type(weight) and "rama_prepro" not in core.scoring.name_from_score_type(weight)):
        lr_terms.set(weight, 0)



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


        self.scorefxn_none = ros().core.scoring.ScoreFunctionFactory().create_score_function("none")
        self.scorefxn_atr = ros().core.scoring.ScoreFunctionFactory().create_score_function("none")
        self.scorefxn_atr.set_weight(ros().core.scoring.fa_atr, 1)
        self.scorefxn_beta = pyro().get_fa_scorefxn()
        self.scorefxn_beta_soft = ros().core.scoring.ScoreFunctionFactory().create_score_function("beta_nov16_soft")


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


def mutations_vs_native(mutant_seq, native_seq):
    mutation_dict = {}
    for seqpos0 in range(len(mutant_seq)):
        if mutant_seq[seqpos0] != native_seq[seqpos0]:
            mutation_dict[seqpos0+1] = mutant_seq[seqpos0]
    return mutation_dict



name3_to_name1 = {
    'ALA':'A',
    'CYS':'C',
    'ASP':'D',
    'GLU':'E',
    'PHE':'F',
    'GLY':'G',
    'HIS':'H',
    'ILE':'I',
    'LYS':'K',
    'LEU':'L',
    'MET':'M',
    'ASN':'N',
    'PRO':'P',
    'GLN':'Q',
    'ARG':'R',
    'SER':'S',
    'THR':'T',
    'VAL':'V',
    'TRP':'W',
    'TYR':'Y'
}

name_used = defaultdict(lambda : 0)

def worst_possible_asp(in_pose, name_no_suffix, out_score_map, out_string_map, suffix):


    monomer_size = in_pose.conformation().chain_end(1)

    in_mutations = mutations_vs_native(in_pose.sequence(), og_pose.sequence())
    in_labels = []
    for seqpos in range(1, in_pose.size()+1):
        for label in in_pose.pdb_info().get_reslabels(seqpos):
            if label.startswith('internal') or label.startswith('backup'):
                in_labels.append((seqpos, label))

    out_stuff = []
    for other_seq, other_labels in zip(all_sequences, all_labels):

        other_mutations = mutations_vs_native(other_seq, og_pose.sequence())

        compatible = True
        for in_seqpos, in_letter in in_mutations.items():
            if in_seqpos not in other_mutations:
                continue
            other_letter = other_mutations[in_seqpos]
            if other_letter == in_letter:
                continue
            compatible = False
            break

        # for seqpos, letter in in_mutations.items():
        #     print(seqpos, letter)
        # for seqpos, letter in other_mutations.items():
        #     print(seqpos, letter)

        if not compatible:
            # print("Not compatible")
            continue

        # Same pdb
        if len(set(list(in_mutations)) | set(list(other_mutations))) == len(in_mutations):
            continue

        pose = og_pose.clone()
        working_seq = list(pose.sequence())

        outer_parts = []
        for iset, (mut_dict, labels) in enumerate([[in_mutations, in_labels], [other_mutations, other_labels]]):

            for position, letter in mut_dict.items():
                # protocols.toolbox.pose_manipulation.repack_this_residue(position, pose, scorefxn, False, letter)
                working_seq[position-1] = letter

            local_tag = []
            for seqpos, label in labels:
                pose.pdb_info().add_reslabel(seqpos, label + str(iset+1))
                local_tag.append(f'{seqpos}{mut_dict[seqpos]}')
            outer_parts.append('_'.join(local_tag))

        rosetta_packer.thread_seq(pose, ''.join(working_seq))

        tag = os.path.basename(args.original_pdb).replace('.pdb', '') + '_' + outer_parts[0] + '-_' + outer_parts[1] + '_v5'

        name_used[tag] += 1
        tag += '_design%i'%name_used[tag]

        print("Generated:", tag)

        score_map = std.map_std_string_double()
        string_map = std.map_std_string_std_string()

        out_stuff.append([pose, tag, score_map, string_map])

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


all_sequences = []
all_labels = []
for pdb in pdbs:
    if ( silent == "" ):
        pose = pose_from_file(pdb)
    else:
        pose = Pose()
        sfd_in.get_structure(pdb).fill_pose(pose)

    all_sequences.append(pose.sequence())
    these_labels = []
    for seqpos in range(1, pose.size()+1):
        for label in pose.pdb_info().get_reslabels(seqpos):
            if label.startswith('internal') or label.startswith('backup'):
                these_labels.append((seqpos, label))
    all_labels.append(these_labels)



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




















