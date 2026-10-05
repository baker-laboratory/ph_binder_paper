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

# import pyRMSD.RMSDCalculator

init("-beta_nov16 -in:file:silent_struct_type binary -keep_input_scores false -mute all"
    " -holes:dalphaball /work/tlinsky/Rosetta/main/source/external/DAlpahBall/DAlphaBall.macgcc"
    )

# setup pH mode
basic.options.set_boolean_option('pH:pH_mode', True)
pose = pose_from_sequence('H')
scorefxn = get_fa_scorefxn()
protocols.toolbox.pose_manipulation.repack_this_residue(1, pose, scorefxn)
basic.options.set_boolean_option('pH:pH_mode', False)




parser = argparse.ArgumentParser()
parser.add_argument("-in:file:silent", type=str, default="")
parser.add_argument("pdbs", type=str, nargs="*")
parser.add_argument("--replicates", type=int, default=1)
parser.add_argument("--low_ph_same_rotamers", action='store_true')
parser.add_argument("--never_pack", action='store_true')


args = parser.parse_args(sys.argv[1:])

pdbs = args.pdbs
silent = args.__getattribute__("in:file:silent")



scorefxn = get_fa_scorefxn()
scorefxn_fa_atr = core.scoring.ScoreFunctionFactory.create_score_function("none")
scorefxn_fa_atr.set_weight(core.scoring.fa_atr, 1)

scorefxn_none = core.scoring.ScoreFunctionFactory.create_score_function("none")
scorefxn_elec = core.scoring.ScoreFunctionFactory.create_score_function("none")
scorefxn_elec.set_weight(core.scoring.fa_elec, 1)


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


scorefxn_hbonds = get_fa_scorefxn()

    
def fix_scorefxn(sfxn, allow_double_bb=False):
    opts = sfxn.energy_method_options()
    opts.hbond_options().decompose_bb_hb_into_pair_energies(True)
    opts.hbond_options().bb_donor_acceptor_check(not allow_double_bb)
    sfxn.set_energy_method_options(opts)

fix_scorefxn(scorefxn_hbonds, allow_double_bb=True)


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

    for seqpos in range(1, pose.size()+1):
        if not his_sub[seqpos]:
            continue
        if do_protonate:
            assert '_P' in pose.residue(seqpos).name(), 'Protonation failed. Do you have -pH_mode True?'
        else:
            assert '_P' not in pose.residue(seqpos).name(), 'Deprotonation failed. Do you have -pH_mode True?'


    basic.options.set_boolean_option('pH:pH_mode', False)

def classify_binder_pos(pose, any_ddg_is_interface=False):
    chainA = core.select.residue_selector.ChainSelector("A")
    chainB = core.select.residue_selector.ChainSelector("B")
    interface_on_A = core.select.residue_selector.NeighborhoodResidueSelector(chainB, 10.0, False)
    interface_on_B = core.select.residue_selector.NeighborhoodResidueSelector(chainA, 10.0, False)
    interface_by_vector = core.select.residue_selector.InterGroupInterfaceByVectorSelector(interface_on_A, interface_on_B)
    interface_by_vector.cb_dist_cut(11)
    interface_by_vector.cb_dist_cut(5.5)
    interface_by_vector.vector_angle_cut(75)
    interface_by_vector.vector_dist_cut(9)

    monomer = pose.split_by_chain()[1]
    sequence = monomer.sequence()
    # dssp = "x" + core.scoring.dssp.Dssp(monomer).get_dssp_secstruct()

    dssp = better_dssp3(monomer)

    sc_neighbors = core.select.util.SelectResiduesByLayer()
    sc_neighbors.use_sidechain_neighbors( True )
    sc_neighbors.compute(pose, "", True)

    atomic_depth = core.scoring.atomic_depth.AtomicDepth( pose, 2.3, False, 0.5 )
    atomic_depth_monomer = core.scoring.atomic_depth.AtomicDepth( monomer, 2.3, False, 0.5 )
    type_set = pose.residue(1).type().atom_type_set()

    probe_size = 2.8
    per_atom_sasa = core.id.AtomID_Map_double_t()
    rsd_sasa = utility.vector1_double()
    core.scoring.calc_per_atom_sasa(pose, per_atom_sasa, rsd_sasa, 2.8, False)

    scorefxn(pose)
    scorefxn(monomer)

    interface_subset = interface_by_vector.apply(pose)

    pose_dats = []

    for seqpos in range(1, monomer.size()+1):

        data = {"ssm_seqpos":str(seqpos)}

        data['sc_neighbors'] = sc_neighbors.rsd_sasa(seqpos)
        data['is_loop'] = dssp[seqpos] == "L"
        data['by_vector'] = interface_subset[seqpos]
        data['dssp'] = dssp[seqpos]
        data['name1'] = monomer.residue(seqpos).name1()

        res = pose.residue(seqpos)
        monomer_res = monomer.residue(seqpos)
        data['depth'] = atomic_depth.calcdepth(res.atom(res.nbr_atom()), type_set)
        data['depth_monomer'] = atomic_depth_monomer.calcdepth(monomer_res.atom(monomer_res.nbr_atom()), type_set)

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
    if any_ddg_is_interface:
        has_ddg = np.abs(pose_df['ddg']) > 0.1
    is_core = ((pose_df['sc_neighbors'] > 5.2) | (pose_df['depth'] - pose_df['depth_monomer'] > 1 )) & (pose_df['sc_sasa'] < 5)

    pose_df['is_interface_core'] = has_ddg & is_core
    pose_df['is_interface_boundary'] = ( has_ddg | pose_df['by_vector'] ) & ~pose_df['is_interface_core']
    pose_df['is_monomer_core'] = ~has_ddg & (pose_df['sc_neighbors'] >= 5.2) & ~pose_df['is_interface_boundary']
    pose_df['is_monomer_boundary'] = ~has_ddg & (pose_df['sc_neighbors'] < 5.2) & (pose_df['sc_neighbors'] >= 2.0) & ~pose_df['by_vector']
    pose_df['is_monomer_surface'] = ~has_ddg & (pose_df['sc_neighbors'] < 2.0) & ~pose_df['by_vector']

    # Make sure all positions are in exactly 1 category
    assert( np.all( pose_df[['is_interface_core', 'is_interface_boundary', 'is_monomer_core', 'is_monomer_boundary',
                            'is_monomer_surface']].sum(axis=1) == 1) )


    return pose_df


def get_simple_is_core(pose):

    sc_neighbors = core.select.util.SelectResiduesByLayer()
    sc_neighbors.use_sidechain_neighbors( True )
    sc_neighbors.compute(pose, "", True)
    is_core = [None]
    for seqpos in range(1, pose.size()+1):
        is_core.append( sc_neighbors.rsd_sasa(seqpos) > 4)

    return is_core


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


def calc_ddg_norepack(pose, scorefxn):
    pose = pose.clone()
    
    close_score = scorefxn(pose)
    ce_close = interface_energy(pose, scorefxn)
    pose = move_chainA_far_away(pose)
    far_score = scorefxn(pose)
    ce_far = interface_energy(pose, scorefxn)

    ddg = close_score - far_score
    ce_ddg = ce_close - ce_far

    #assert abs(ddg - ce_ddg) < 0.1

    return ddg



def interface_energy(pose, scorefxn, charge_only=False, charge_his_only=False):

    weights = scorefxn.weights()
    gr = pose.energies().energy_graph()
    monomer_size = pose.conformation().chain_end(1)
    dg = 0
    for seqpos in range(1, monomer_size+1):
        name1 = pose.residue(seqpos).name1()
        for seqpos2 in range(monomer_size+1, pose.size()+1):
            name2 = pose.residue(seqpos2).name1()

            if charge_only:
                if name1 not in 'DERKH':
                    continue
                if name2 not in 'DERKH':
                    continue
            if charge_his_only:
                good = False
                if name1 in 'DERKH' and name2 == 'H':
                    good = True
                if name2 in 'DERKH' and name1 == 'H':
                    good = True
                if not good:
                    continue


            edge = gr.find_edge(seqpos, seqpos2)
            if not edge:
                continue
            dg += edge.dot(weights)

    return dg


def calc_ddg_norepack_charge_only(pose, scorefxn, charge_only=False, charge_his_only=False):
    pose = pose.clone()
    
    scorefxn(pose)
    close_score = interface_energy(pose, scorefxn, charge_only=charge_only, charge_his_only=charge_his_only )
    pose = move_chainA_far_away(pose)
    scorefxn(pose)
    far_score = interface_energy(pose, scorefxn, charge_only=charge_only, charge_his_only=charge_his_only )

    return close_score - far_score


def mean_reject_edges(array):
    if len(array) < 3:
        return np.mean(array)
    else:
        return np.mean(list(sorted(list(array)))[1:-1])



def rotamer_charged_atoms(rotamer):
    name1 = rotamer.name1()
    atoms = None
    if ( name1 == "D" ):
        atoms = ['OD1', 'OD2']
    if ( name1 == "E" ):
        atoms = ['OE1', 'OE2']
    if ( name1 == "K" ):
        atoms = ['NZ']
    if ( name1 == "R" ):
        atoms = ['NE', 'NH1', 'NH2']
    if ( name1 == "H" ):
        atoms = [rotamer.atom_name(x).strip() for x in range(rotamer.first_sidechain_atom(), rotamer.nheavyatoms()+1) 
                                                                            if rotamer.atom_name(x).strip().startswith('N')]
    assert( not atoms is None)

    xyzs = np.zeros((len(atoms), 3))
    for i in range(len(atoms)):
        xyzs[i] = from_vector(rotamer.xyz(atoms[i]))
    return xyzs

def from_vector(xyz_vector):
    return np.array([xyz_vector.x, xyz_vector.y, xyz_vector.z])

def get_CBs(pose):
    CBs = []
    for seqpos in range(1, pose.size()+1):
        CBs.append(from_vector(pose.residue(seqpos).nbr_atom_xyz()))
    CBs = np.array(CBs)
    return CBs


def get_salt_score(pose, salt_weight=-2, low_pH=False):

    CBs = get_CBs(pose)

    monomer_size = pose.conformation().chain_end(1)
    salt_ddg = 0
    for seqpos in range(1, monomer_size+1):
        salt_ddg += salt_weight * get_salt_score_ind(pose, seqpos, CBs, low_pH=low_pH)

    return salt_ddg

def get_salt_score_ind(pose, seqpos, CBs, low_pH=False):
    sequence = pose.sequence()

    monomer_size = pose.conformation().chain_end(1)
    is_target = np.zeros(pose.size(), bool)
    is_target[monomer_size:] = True

    main_letter = pose.residue(seqpos).name1()
    charge_letters = 'DERKH' if low_pH else 'DERK'
    if ( main_letter not in charge_letters ):
        return 0

    we_are_negative = main_letter in "DE"

    cb_dist_from_us = np.linalg.norm( CBs - CBs[seqpos-1], axis=-1 )
    close_enough_to_count = cb_dist_from_us < 16

    if ( we_are_negative ):
        ok_letter = (np.array(list(sequence)) == "K") | (np.array(list(sequence)) == "R") | (np.array(list(sequence)) == "H")
        is_potential_partner = is_target & close_enough_to_count & ok_letter

        anti_letter = (np.array(list(sequence)) == "D") | (np.array(list(sequence)) == "E")
        is_anti_partner = is_target & close_enough_to_count & anti_letter
    else:
        ok_letter = (np.array(list(sequence)) == "D") | (np.array(list(sequence)) == "E")
        is_potential_partner = is_target & close_enough_to_count & ok_letter

        anti_letter = (np.array(list(sequence)) == "R") | (np.array(list(sequence)) == "K") | (np.array(list(sequence)) == "H")
        is_anti_partner = is_target & close_enough_to_count & anti_letter


    our_atoms = rotamer_charged_atoms(pose.residue(seqpos))

    partners = np.concatenate((np.where(is_potential_partner)[0]+1, -(np.where(is_anti_partner)[0]+1)))

    salt_before_weight = 0

    for partner in partners:
        multiplier = np.sign(partner)
        partner = abs(partner)

        their_atoms = rotamer_charged_atoms(pose.residue(partner))

        closest = np.min(np.linalg.norm( our_atoms[:,None] - their_atoms[None,:], axis=-1 ))

        score = 0
        if ( closest < 7 ):
            score = 1
        elif ( closest < 9 ):
            score = 1 - (closest - 7) / (9 - 7)
        else:
            score = 0

        salt_before_weight += score * multiplier

    # if ( seqpos == 102 and main_letter == "R" ):
    #     import IPython
    #     IPython.embed()

    # print("Salt before weight %.2f"%(salt_before_weight))
    return salt_before_weight





the_locals = None

def worst_possible_asp(pose, name_no_suffix, out_score_map_real, out_string_map, suffix):


    monomer_size = pose.conformation().chain_end(1)

    is_core = get_simple_is_core(pose)

    score_maps = []
    for replicate in range(args.replicates):

        out_score_map = {}
        ddgs = []
        salt_parts = []
        ddgs_w_salt = []
        ddg_elecs = []
        ddg_elecs_charged = []
        ddg_elecs_chargedH = []
        net_charges = []

        for low_pH in [False, True]:

            pre_packed = pose.clone()
            protonate_histidines(pose, low_pH)

            if args.low_ph_same_rotamers and low_pH:
                for seqpos in range(1, pose.size()):
                    for chi in range(1, pose.residue(seqpos).nchi()+1):
                        pose.set_chi(chi, seqpos, pre_packed.residue(seqpos).chi(chi))
                for seqpos in range(1, pose.size()):
                    for chi in range(1, pose.residue(seqpos).nchi()+1):
                        assert np.isclose(pose.residue(seqpos).chi(chi), pre_packed.residue(seqpos).chi(chi))
            else:
                if not args.never_pack:
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

            for cat in cats2:
                out_score_map[prefix + cat] = 0

            for seqpos in range(1, pose.size()+1):
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


                print(f'{prefix} {seqpos:3d} {pose.pdb_info().chain(seqpos)} {cat:7s} cross hbonds: {cross_hbonds}')

                assert prefix + cat in out_score_map
                out_score_map[prefix + cat] += cross_hbonds


            monomer_seq = pose.sequence()[:monomer_size]
            net_charge = monomer_seq.count('R') + monomer_seq.count('K') - monomer_seq.count('D') - monomer_seq.count('E')
            if low_pH:
                net_charge += monomer_seq.count('H')

            ddg_pose = pose.clone()
            rosetta_packer.beta_pack(ddg_pose)
            ddgs.append(calc_ddg_norepack(ddg_pose, scorefxn))
            salt_parts.append(get_salt_score(ddg_pose))
            ddgs_w_salt.append(ddgs[-1] + salt_parts[-1])
            ddg_elecs.append(calc_ddg_norepack(ddg_pose, scorefxn_elec))
            ddg_elecs_charged.append(calc_ddg_norepack_charge_only(ddg_pose, scorefxn_elec, charge_only=True))
            ddg_elecs_chargedH.append(calc_ddg_norepack_charge_only(ddg_pose, scorefxn_elec, charge_his_only=True))
            net_charges.append(net_charge)

            # pose.dump_pdb(f"{low_pH}_pH.pdb")


        ph_score = calc_ph_score(out_score_map)
        out_score_map['ph_score'] = ph_score

        out_score_map['ddg_soft_lowph_vs_high'] = ddgs[1] - ddgs[0]
        out_score_map['ddg_soft_w_salt_lowph_vs_high'] = ddgs_w_salt[1] - ddgs_w_salt[0]
        out_score_map['salt_part_lowph_vs_high'] = salt_parts[1] - salt_parts[0]
        out_score_map['ddg_elec_lowph_vs_high'] = ddg_elecs[1] - ddg_elecs[0]
        out_score_map['ddg_elec_lowph_vs_high_charged_only'] = ddg_elecs_charged[1] - ddg_elecs_charged[0]
        out_score_map['ddg_elec_lowph_vs_high_charged_to_HIS_only'] = ddg_elecs_chargedH[1] - ddg_elecs_chargedH[0]
        out_score_map['net_charge_lowph_vs_high'] = net_charges[1] - net_charges[0]

        out_score_map['ddg_soft_highph'] = ddgs[0]
        out_score_map['ddg_soft_lowph'] = ddgs[1]
        out_score_map['ddg_soft_w_salt_highph'] = ddgs_w_salt[0]
        out_score_map['ddg_soft_w_salt_lowph'] = ddgs_w_salt[1]
        out_score_map['salt_part_highph'] = salt_parts[0]
        out_score_map['salt_part_lowph'] = salt_parts[1]
        out_score_map['net_charge_highph'] = net_charges[0]
        out_score_map['net_charge_lowph'] = net_charges[1]


        score_maps.append(out_score_map)

    for key in score_maps[0]:
        out_score_map_real[key] = mean_reject_edges([d[key] for d in score_maps])

    print('pH SCORE:', out_score_map_real['ph_score'])
    print('ddg_soft_lowph_vs_high: %5.1f ddg_elec_lowph_vs_high: %5.1f ddg_elec_lowph_vs_high_charged_only: %5.1f ddg_elec_lowph_vs_high_charged_to_HIS_only: %5.1f'%(out_score_map_real['ddg_soft_lowph_vs_high'], out_score_map_real['ddg_elec_lowph_vs_high'],
        out_score_map_real['ddg_elec_lowph_vs_high_charged_only'], out_score_map_real['ddg_elec_lowph_vs_high_charged_to_HIS_only']))



    return None






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


        core.io.raw_data.ScoreMap.add_arbitrary_score_data_from_pose( pose, score_map)
        core.io.raw_data.ScoreMap.add_arbitrary_string_data_from_pose( pose, string_map)

        sfd.write_pose( pose, score_map, name_no_suffix, string_map)
        if (out_pose != None):

            if ( silent == "" ):
                out_pose.dump_pdb(name_no_suffix + ".pdb")
            else:
                silent_name = "out.silent"
                sfd_out = core.io.silent.SilentFileData( silent_name, False, False, "binary", core.io.silent.SilentFileOptions())
                struct = sfd_out.create_SilentStructOP()
                struct.fill_struct(out_pose, name_no_suffix)
                sfd_out.add_structure(struct)
                sfd_out.write_all(silent_name, False)


        seconds = int(time.time() - t0)

        print("protocols.jd2.JobDistributor: " + name_no_suffix + " reported success in %i seconds"%seconds)

    # except Exception as e:
    #     print("Error!!!")
    #     print(e)



# if ( silent != "" ):
#     sfd_out.write_all("out.silent", False)




















