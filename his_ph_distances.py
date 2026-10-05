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




parser = argparse.ArgumentParser()
parser.add_argument("-in:file:silent", type=str, default="")
parser.add_argument("pdbs", type=str, nargs="*")
parser.add_argument("--num_to_look_at", type=int, default=3)
parser.add_argument("--write_other_seqpos", action='store_true')
parser.add_argument("--use_CB", action='store_true')


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




the_locals = None

def worst_possible_asp(pose, name_no_suffix, out_score_map, out_string_map, suffix):


    monomer_size = pose.conformation().chain_end(1)

    indices = []
    atoms = []
    name1s = []

    for seqpos in range(1, pose.size()):
        res = pose.residue(seqpos)
        if res.name1() not in 'DERKH':
            continue
        indices.append(seqpos)
        name1s.append(res.name1())
        my_atoms = []
        if args.use_CB:
            my_atoms.append(nup.from_vector(res.nbr_atom_xyz()))
        else:

            for iatom in range(5, res.nheavyatoms()+1):
                if res.atom_type(iatom).element() != 'C':
                    my_atoms.append(nup.from_vector(res.xyz(iatom)))
        atoms.append(np.array(my_atoms))

    indices = np.array(indices)
    name1s = np.array(name1s)

    closest_approach = np.zeros((len(indices), len(indices)))

    for ires1 in range(len(indices)):
        for ires2 in range(len(indices)):
            closest = np.linalg.norm(atoms[ires1][:,None] - atoms[ires2][None,:], axis=-1).min()
            closest_approach[ires1, ires2] = closest


    close_enough = closest_approach < 10

    is_monomer = indices <= monomer_size
    is_his = name1s == 'H'
    is_neg = (name1s == 'E') | (name1s == 'D')
    is_pos = (name1s == 'K') | (name1s == 'R')

    cross_chain = is_monomer[:,None] != is_monomer[None,:]
    his_neg = (is_his[:,None] & is_neg[None,:]) & cross_chain & close_enough
    his_pos = (is_his[:,None] & is_pos[None,:]) & cross_chain & close_enough
    his_his = (is_his[:,None] & is_his[None,:]) & cross_chain & close_enough


    his_interesting = his_neg | his_pos | his_his

    seqpos_is_interesting = his_interesting.any(axis=-1)
    seqpos_closest_interesting = np.full((len(indices),), 1000, dtype=float)
    wh_interesting = np.where(his_interesting)
    for ihis, ielse in zip(*wh_interesting):
        seqpos_closest_interesting[ihis] = np.min([seqpos_closest_interesting[ihis], closest_approach[ihis, ielse]])

    seqpos_monomer_interesting = seqpos_is_interesting & is_monomer
    seqpos_target_interesting = seqpos_is_interesting & ~is_monomer



    for on_binder in [True, False]:
        prefix = 'binder' if on_binder else 'target'
        our_interesting = seqpos_monomer_interesting if on_binder else seqpos_target_interesting
        wh_local_interesting = np.where(our_interesting)[0]
        our_closests = seqpos_closest_interesting[our_interesting]
        our_argsort = np.argsort(our_closests)

        sorted_wh = wh_local_interesting[our_argsort]

        for iour in range(args.num_to_look_at):

            if iour >= len(sorted_wh):
                ilookup = -1
                our_seqpos = np.nan
            else:
                ilookup = sorted_wh[iour]
                our_seqpos = indices[ilookup]

            out_score_map[f'{prefix}_his{iour+1}_seqpos'] = our_seqpos

            for label, mask in zip(['his', 'neg', 'pos'], [his_his, his_neg, his_pos]):
                our_mask = mask[ilookup]
                wh_our_mask = np.where(our_mask)[0]
                dists_our_mask = closest_approach[ilookup][our_mask]
                argsort = np.argsort(dists_our_mask)

                wh_our_mask = wh_our_mask[argsort]
                dists_our_mask = dists_our_mask[argsort]

                for iother in range(args.num_to_look_at):
                    if iother < len(wh_our_mask) and ilookup != -1:
                        seqpos_other = indices[wh_our_mask[iother]]
                        dist_other = dists_our_mask[iother]
                        assert pose.residue(our_seqpos).name1() == 'H'
                        if label == 'his':
                            assert pose.residue(seqpos_other).name1() == 'H'
                        elif label == 'pos':
                            assert pose.residue(seqpos_other).name1() in 'RK'
                        else:
                            assert pose.residue(seqpos_other).name1() in 'DE'
                        print('%s -- %3iH -- %5.1f to %3i%s -- %s'%(prefix, our_seqpos, dist_other, seqpos_other, pose.residue(seqpos_other).name1(), label))
                    else:
                        seqpos_other = np.nan
                        dist_other = np.nan
                    out_score_map[f'{prefix}_his{iour+1}_{label}{iother+1}_dist'] = dist_other
                    if args.write_other_seqpos:
                        out_score_map[f'{prefix}_his{iour+1}_{label}{iother+1}_seqpos'] = seqpos_other

                





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




















