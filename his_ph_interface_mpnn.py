#!/usr/bin/env python
# Calls stock ProteinMPNN (https://github.com/dauparas/ProteinMPNN) as an external program for his_ph_interface_design.py
#
# ProteinMPNN is found in this order:
#   1. the PROTEIN_MPNN_PATH environment variable: folder containing protein_mpnn_run.py (a clone of the ProteinMPNN repo)
#   2. a ProteinMPNN folder next to this script (e.g. a git submodule)
# Other environment variables:
#   PROTEIN_MPNN_PYTHON  -- (optional) python that has torch installed. Defaults to the python running this script.
#                           This lets ProteinMPNN live in a different environment than pyrosetta.

import os
import sys
import json
import shutil
import subprocess
import tempfile

import numpy as np


def find_mpnn_path():
    '''Returns the folder with protein_mpnn_run.py or None'''
    candidates = [os.environ.get('PROTEIN_MPNN_PATH'), os.path.join(os.path.dirname(os.path.abspath(__file__)), 'ProteinMPNN')]
    for path in candidates:
        if path and os.path.exists(os.path.join(path, 'protein_mpnn_run.py')):
            return path
    return None


def run_stock_mpnn(design_chains, pdb_str, best_of_n=1, return_all=True, sampling_temps=(0.1,), bias_by_res_d=None,
                    fixed_pos_list_by_chain_d=None, tied_positions_list=None):
    '''
    Returns (sequences, scores). One entry for every temperature and for every best_of_n.
    Sequences contain the designed chains only, in chain order.

    design_chains -- 'A' or 'A C' (space separated)
    bias_by_res_d -- {chain: array of shape (length, 21)}, alphabet is ACDEFGHIKLMNPQRSTVWYX
    fixed_pos_list_by_chain_d -- {chain: [1-indexed positions]}
    tied_positions_list -- [{'A':[1], 'C':[1]}, ...]
    '''

    mpnn_path = find_mpnn_path()
    assert mpnn_path is not None, "Can't find ProteinMPNN. Set PROTEIN_MPNN_PATH to a clone of https://github.com/dauparas/ProteinMPNN"
    mpnn_python = os.environ.get('PROTEIN_MPNN_PYTHON', sys.executable)
    run_script = os.path.join(mpnn_path, 'protein_mpnn_run.py')

    name = 'input'
    tmpdir = tempfile.mkdtemp(prefix='mpnn_')
    try:
        pdb_path = os.path.join(tmpdir, name + '.pdb')
        with open(pdb_path, 'w') as f:
            f.write(pdb_str)

        cmd = [mpnn_python, run_script,
                '--pdb_path', pdb_path,
                '--pdb_path_chains', design_chains,
                '--out_folder', tmpdir,
                '--num_seq_per_target', str(best_of_n),
                '--batch_size', '1',
                '--sampling_temp', ' '.join(str(x) for x in sampling_temps),
                '--omit_AAs', 'C',
                '--suppress_print', '1',
                ]

        def add_jsonl(flag, filename, data):
            path = os.path.join(tmpdir, filename)
            with open(path, 'w') as f:
                f.write(json.dumps({name: data}) + '\n')
            cmd.extend([flag, path])

        if bias_by_res_d is not None:
            add_jsonl('--bias_by_res_jsonl', 'bias_by_res.jsonl', {chain: np.asarray(bias).tolist() for chain, bias in bias_by_res_d.items()})
        if fixed_pos_list_by_chain_d is not None:
            add_jsonl('--fixed_positions_jsonl', 'fixed.jsonl', {chain: [int(x) for x in positions] for chain, positions in fixed_pos_list_by_chain_d.items()})
        if tied_positions_list is not None:
            add_jsonl('--tied_positions_jsonl', 'tied.jsonl', tied_positions_list)

        # One thread is fastest (and keeps torch from grabbing every core)
        env = dict(os.environ)
        env.setdefault('OMP_NUM_THREADS', '1')
        result = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, universal_newlines=True, env=env)
        if result.returncode != 0:
            raise RuntimeError("ProteinMPNN failed:\n" + ' '.join(cmd) + "\n" + result.stdout + "\n" + result.stderr)

        with open(os.path.join(tmpdir, 'seqs', name + '.fa')) as f:
            lines = [x.strip() for x in f if x.strip()]
    finally:
        shutil.rmtree(tmpdir, ignore_errors=True)

    # The first record is the input sequence, then the samples
    seqs = []
    scores = []
    for header, seq in zip(lines[2::2], lines[3::2]):
        fields = dict(x.strip().split('=', 1) for x in header.lstrip('>').split(',') if '=' in x)
        seqs.append(seq.replace('/', ''))
        scores.append(float(fields['score']))

    assert len(seqs) == best_of_n * len(sampling_temps), (len(seqs), best_of_n, sampling_temps)

    return seqs, scores
