# pH-dependent protein binders with histidines

Scripts for turning existing protein binders into pH-dependent binders by placing histidines (HIS).

From our paper, we described two separate ways to create pH dependence:

* One where the binder structure is disrupted at low pH. We internally termed this the "his_ph_exploder" pipeline, and it can be found in the [exploder pipeline](#exploder-pipeline) section.
* The other is described as "his_ph_interface" ([interface pipeline](#interface-pipeline)) and works by creating cross-interface features that are destabilized at low pH.

The two pipelines are separate but may be combined on the same design if you want.

## Installation

```bash
git clone --recursive git@github.com:baker-laboratory/ph_binder_paper.git   # --recursive pulls in ProteinMPNN
conda create -n his_ph python=3.12
conda activate his_ph
pip install npose scipy numpy numba pandas biopython torch
# plus PyRosetta, see https://www.pyrosetta.org
```

See `environment.txt` for details. ProteinMPNN (`ProteinMPNN/`, a git submodule pinned to the tested commit) needs torch. If you cloned without `--recursive`, run `git submodule update --init`. Set `OMP_NUM_THREADS=1` when running ProteinMPNN, it is fastest.

An oracle such as AF2 or AF3 is used by both pipelines and is not included.

### Usage conventions

Every script takes pdbs or a silent file:

```bash
python script.py *.pdb
python script.py -in:file:silent input.silent
```

Outputs are PDBs (or `out.silent`) plus a `score.sc` scorefile. Run any script with `-h` for its options.

## Exploder pipeline

The general idea is to put h-bonding histidines in the core of pre-existing binders so that they will become pH dependent.

**Two regimes: HIS-HIS h-bonds vs. HIS next to a positive charge.** While the original hypothesis was that pairs of histidines h-bonding to each other would provide the strongest effect through steric clashing of the hydrogens, it was later found that the positive-positive charge effect is perhaps the stronger force. To that end, these scripts come with options to not worry about making h-bonds but rather seek to put a HIS next to an ARG or a LYS (`--other_aas RK --dist_cutoff_mode <Angstroms>`, usually with `--interface_mode`). If one were doing a screening campaign, it may be wise to order designs under both regimes.

This protocol was originally run on several combo variants identified from yeast display, so there were multiple copies of the same binder. Multiple known-good copies of the same binder is a good starting point, but without them you can follow everything except the threading step (6a).

```
input binder
  -> 1.  add HIS pairs          his_ph_exploder_v3.py
  -> 2.  add backup mutations   his_ph_exploder_backup_v3.py
  -> 3.  design the region      ProteinMPNN + his_ph_exploder_design_v2.py
  -> 4.  oracle (AF2 / AF3)
  -> 5.  check the hbnets       his_ph_exploder_evaluate.py
  -> 6a. (optional) thread known-good networks   his_ph_exploder_thread_mutations.py
  -> 6b. (optional) merge two networks           his_ph_exploder_merge_mutations.py
  -> 7.  choose designs to order
```

### 1. Add the first histidine pairs

```bash
python his_ph_exploder_v3.py *.pdb
```

This script only really cares about your backbone. If you don't get any outputs, it may not be possible to make this binder pH dependent in this way. If you really want this to work:

* If your binder is experimentally validated, try predicting the structure with as many oracles as you can find (including Rosetta FastRelax).
* If your binder is computational, perhaps try some partial diffusion.

The script also reports the solvent accessibility of your HIS. The most important column is `hbond_depth`, which tells you how deep into the protein your HIS is. Depth < 4 is on the surface and you should probably throw those away. Energy worse than -0.5 is going to get removed at a later step, however it is possible it will get better after design. `max_sc_sasa` is the SASA of the original sidechain, so it is less useful than you'd hope.

Flags:

* `--other_aas RK` also tries HIS-ARG and HIS-LYS pairs.
* `--interface_mode` only looks for pairs across the interface.
* `--two_sided_design` (use with `--interface_mode`) allows the target to mutate.
* `--dist_cutoff_mode <Angstroms>` is for HIS next to a positive charge. It does not look for h-bonds. Instead, two non-carbon atoms within this many Angstroms of each other count as a pair. Likely you will also use `--other_aas RK` and `--interface_mode` with this flag. The score column is then `min_distance` instead of `energy`.
* `--only_allow_positions 1,2,3` restricts the search to these (1-indexed) positions.

### 2. Add backup mutations

Run on the output of step 1:

```bash
python his_ph_exploder_backup_v3.py *.pdb
```

This looks for positions near your HIS pairs where another HIS or a SER could be placed to increase the size of the h-bond network. Whether this helps or hurts experimentally is unknown, but it generates more examples to test.

* `--backup_aas HSRK` additionally backs up with R and K (default `HS`).
* `--allow_interface` also allows backup positions on the other side of the interface.
* `--input_doesnt_hbond` does not require the input HIS pair to h-bond each other. The backup search then starts from an empty atom set. By default the input pair must h-bond.
* `--two_sided_design` allows the target to mutate.

### 3. Design the region around the mutations

Collect the outputs of **both steps 1 and 2**.

First, generate the ProteinMPNN conditional probabilities for each pdb:

```bash
export OMP_NUM_THREADS=1
python ProteinMPNN/protein_mpnn_run.py --pdb_path_chains A --conditional_probs_only 1 \
    --out_folder ./mpnn_probs --pdb_path pdb_of_interest.pdb
```

Put the full path of all the `.npz` files into a single list file:

```bash
readlink -f mpnn_probs/conditional_probs_only/*.npz > mpnn_probs.list
```

Then run the design script on all of the outputs of steps 1 and 2:

```bash
python his_ph_exploder_design_v2.py -in:file:silent step1and2outputs.silent --list_of_mpnn_prob_npz mpnn_probs.list
```

This script is slow, up to about 40 minutes per input, so plan for that when submitting. `--modes` selects what it does (default `regular,force_touching`), and `--two_sided_design` allows the target to mutate.

### 4. Oracle your designs

Run the outputs of step 3 through your favorite oracle to get new structures. If you are working with a validated binder, make sure that the `pae_interaction` does not get worse by more than about 3 and that the interface RMSD to the original is not worse than about 2 A. (`pae_interaction` comes from your oracle. `his_ph_exploder_energy_evaluate.py` reports an `interface_rmsd`, see step 7.)

It is critical that the oracle maintains the PDBInfoLabels in your pdbs. If your oracle does not support this, they are just short lines of text at the bottom of the pdb that you can copy over. Split your pdbs by length as usual.

### 5. Check that the hbnets still look good

Passing the oracle does not mean the HIS are still h-bonding to each other:

```bash
python his_ph_exploder_evaluate.py *.pdb
```

Columns:

* `his_hbond`: energy of the primary HIS h-bond (out of -2.1)
* `mean_backup_e`: energy of the backup h-bond, if it exists
* `n_backups`: how many backup residues there are
* `n_satisfied_backups`: how many backup residues are satisfied

In general, call a design a success if `his_hbond < -0.5` and `n_backups == n_satisfied_backups`. You could order these at this point.

### 6a. (Optional) Thread known-good networks onto similar designs

If you have several binders with the same backbone (for example, combo variants from yeast display), you can plant networks that you know are good into the other backbones to get more outputs. Steps 1 and 2 should in principle find these networks anyway, but this can add some.

Make a list of the networks, one per line. The first two entries are the primary HIS and the third, optional, entry is the backup residue:

```
11H,38H
15H,91H
15H,91H 79S
45H,65H
45H,65H 8H
```

(Easy to produce with `sed` on the names of your successful outputs.) Then, on your original inputs:

```bash
python his_ph_exploder_thread_mutations.py *.pdb --mutation_thread_list networks.list
```

The outputs go straight into step 3.

### 6b. (Optional) Mix and match networks

This creates new designs containing two networks from your good designs (up to 6 HIS in a single design). The success rate is low, but these designs are presumably the most likely to be pH responsive.

The input is a little different than usual: group your inputs by the original parent pdb (the binder with no HIS), then run:

```bash
python his_ph_exploder_merge_mutations.py good_designs_of_this_parent*.pdb --original_pdb parent.pdb
```

The outputs go into step 3. At step 5, add `--num_sets 2` and make sure that both set 1 and set 2 pass the criteria in the resulting scorefile. Most will fail, but a few should succeed.

### 7. Choosing designs to order

What makes a good design:

* The histidines make a good h-bond
* The histidines are buried
* `pae_interaction` is close to the original design
* The interface RMSD to the original design is < 1
* The histidines are on separate secondary structure elements near the interface
* The h-bonds remain after FastRelax

`his_ph_exploder_energy_evaluate.py` relaxes each design at low and high pH and reports the energy difference (`delta_score_ph`), `interface_rmsd` and `monomer_rmsd`. It is slow (minutes per input) and the energy difference is noisy, so use it with caution. At the least, it can be used to FastRelax your designs, after which you can recheck the hbnets with step 5. Use `--nstruct` to set the number of relaxations.

### Visualizing in PyMOL

The names of the designs tell you where the HIS are. This PyMOL snippet gives a good visualization:

```
pdbdircyclerlite
set_cycler_command "stix; cmd.select('sele', 'resi ' + '+'.join([x[:-1] for x in cmd.get_object_list()[0].replace('-', '').split('_') if x[0].isnumeric() and x[-1].isalpha()])); color yellow, sele; zoom sele, 10;  deselect; util.cba('all')"
```

## Interface pipeline

This pipeline designs HIS **at the interface** so that they h-bond the target at high pH and lose that h-bond (and pick up positive charge) at low pH.

```
binder:target complex (no HIS needed)
  -> 1. his_ph_interface_design.py   picks HIS positions, designs with ProteinMPNN, scores pH sensitivity
  -> 2. oracle (AF2 / AF3)
  -> 3. his_ph_interface_hbonds.py   re-scores the predicted structures
```

The complex must have the binder as chain A and the target as chain B. The binder does not need to contain any HIS.

### 1. Design

```bash
python his_ph_interface_design.py binder_on_target.pdb
```

What the script does:

1. Finds binder positions that could take a HIS that accepts an h-bond from the target. Each candidate gets an initial score: 6 for a core backbone acceptor, 3 for a core sidechain acceptor, 2 for a surface backbone acceptor and 1 for a surface sidechain acceptor.
2. First round: draws `-first_round_n_sets` sets of `-his_in_groups_of` positions (weighted by that score), places HIS there, designs the rest of the binder with ProteinMPNN, repacks in Rosetta and gives each set a pH score.
3. Second round: draws sets again from the positions that scored well (`-second_round_n_sets`).
4. Keeps designs with `ph_score > min_his_score`, at most `-num_per_input` per input.

The pH score rewards HIS that make cross-interface h-bonds at high pH (neutral HIS, more for core positions and for backbone partners) and penalizes HIS that still h-bond the target at low pH (protonated HIS). Higher is better.

Outputs are `out.silent` and `score.sc` (columns `ph_score`, `sap_score`, `ss_schain_atr`, `ddg_norepack_soft`, `contact_patch`). Designs are named `<input>_<positions>H_..._mpnn<nnn>`, for example `Tie2_tw4104_11H_15H_26H_50H_mpnn003`, and each kept HIS gets a `ph_score:<x>` PDBInfoLabel. If no good HIS are found for an input the script reports "There were no outputs?" and moves on.

Main options:

| Option | Default | Meaning |
|---|---|---|
| `-his_in_groups_of` | 3 | HIS per trial set |
| `-first_round_n_sets`, `-second_round_n_sets` | 10, 10 | Number of sets tried in each round |
| `-min_his_score` | 5 | Minimum pH score to output |
| `-num_per_input` | 20 | Maximum outputs per input |
| `-mpnn_seqs`, `-mpnn_temps` | 1, `"0.001 0.01 0.1"` | ProteinMPNN sampling |
| `-hbnet_lock_identities`, `-hbnet_lock_identities_strict` | off | Lock the identities of residues in the existing h-bond network |

Run `python his_ph_interface_design.py -h` for the rest.

ProteinMPNN is run as an external program. It is found automatically in the `ProteinMPNN/` submodule. To use another copy set `PROTEIN_MPNN_PATH`, and if that copy lives in a different python environment (the one with torch), set `PROTEIN_MPNN_PYTHON`. It runs with `OMP_NUM_THREADS=1` unless you set the variable yourself.

### 2. Oracle

Predict the designed complexes with AF2, AF3 or similar (split inputs by length as usual). Compare `pae_interaction` and the interface RMSD against the original binder:target complex and drop designs that moved.

### 3. Re-score the predicted structures

```bash
python his_ph_interface_hbonds.py oracle_outputs/*.pdb [--replicates 3] [--low_ph_same_rotamers] [--never_pack]
```

For each complex this protonates all HIS (low pH) and deprotonates them (high pH), repacks each state and writes one row per input to `score.sc`. Columns:

* `ph_score`: the same score as in the design step, now on the predicted structure. This is the number to filter on.
* `high_pH_*`, `low_pH_*`: counts of cross-interface HIS h-bonds by category (`core_`/`surf_`; `A` acceptor, `D` donor; `_bb` when the partner is backbone)
* `ddg_soft_highph`, `ddg_soft_lowph`, `ddg_soft_lowph_vs_high`: soft binding energy at each pH. A positive `lowph_vs_high` means binding is weaker at low pH, which is what you want.
* `ddg_elec_lowph_vs_high` (plus `_charged_only` and `_charged_to_HIS_only`), `ddg_soft_w_salt_*`, `salt_part_*`: electrostatic and salt bridge terms
* `net_charge_highph`, `net_charge_lowph`

`--replicates N` repeats the packing and averages the replicates, `--low_ph_same_rotamers` reuses the high pH rotamers at low pH, and `--never_pack` skips repacking. The script does not write pdbs.

This scorer does not depend on the design script, so it can also be run on the original complex for a baseline of the `ddg_*` columns (`ph_score` is 0 when the binder has no HIS) or on any other complex that contains HIS.

### Other scripts

Standalone screening scripts for the same kind of input:

* `his_ph_interface_hbond_finder.py` lists binder positions that accept an h-bond from the target and classifies them as `core_bb_acc`, `surf_bb_acc`, `core_sc_acc` or `surf_sc_acc` (position lists joined by `+` in the scorefile). Useful to see where HIS could go without running the design.
* `his_ph_interface_distances.py` finds, for each binder HIS, the closest approach to target HIS, ARG/LYS (`pos`) and ASP/GLU (`neg`) residues within 10 A, using non-carbon sidechain atoms (or CB with `--use_CB`). Columns look like `binder_his1_pos1_dist` and `binder_his1_pos1_seqpos`. `--num_to_look_at` (default 3) sets how many are reported and `--write_other_seqpos` adds positions. This is the "HIS near a positive charge" idea as a screen rather than a design script.

## Script index

| Script | Pipeline | Step |
|---|---|---|
| `his_ph_exploder_v3.py` | exploder | 1 |
| `his_ph_exploder_backup_v3.py` | exploder | 2 |
| `his_ph_exploder_design_v2.py` | exploder | 3 |
| `his_ph_exploder_evaluate.py` | exploder | 5 |
| `his_ph_exploder_thread_mutations.py` | exploder | 6a |
| `his_ph_exploder_merge_mutations.py` | exploder | 6b |
| `his_ph_exploder_energy_evaluate.py` | exploder | 7 |
| `his_ph_interface_design.py` | interface | 1 |
| `his_ph_interface_hbonds.py` | interface | 3 |
| `his_ph_interface_hbond_finder.py` | interface | screening |
| `his_ph_interface_distances.py` | interface | screening |
| `his_ph_interface_mpnn.py` | interface | helper that runs ProteinMPNN |

## License

MIT, see `LICENSE`.
