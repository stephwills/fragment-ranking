## fragment-ranking

A method for ranking fragments by how much novel information they give about protein targets in fragment screens. When using the results of fragment screens on many diverse targets, this method has been shown to select a set of functionally diverse fragments that can get information more efficiently from new targets. 

# Instructions for use

Installation:

1. clone repository
2. :construction:


Required input files:
1. Fragment library, in sdf format.
2. Structures of targets with fragments bound, in pdb format. Structures should not have any atoms or residues missing.
3. JSON dictionary of pdb file with associated fragment (in SMILES string), for example : {"mArh-x1018.pdb": "Clc1ccc2nnnn2n1",...}

Test data is provided.

Run:
1. Generate interaction fingerprints from structures. Either residue-level or atomic-level interaction fingerprints can be used. 
```bash
$ python src/generate_IFPs.py -IFP [atomic/residue] -sdf [path_to_fragment_library.sdf] -exps [path_to_experiments.json] -pdbs [path_to_pdb_files]
```
  for example:
```bash
$ python src/generate_IFPs.py -IFP atomic -sdf data/library.sdf -exps data/experiments_mArh.json -pdbs data/structures/
```
2. Rank fragments
```bash
$ python src/rank_fragments.py -sdf [path_to_fragment_library.sdf] -bits data/smiles_bits_[atomic/residue].json -o [output_file]
```
for example:
```bash
$ python src/rank_fragments.py -sdf data/library.sdf -bits data/smiles_bits_atomic.json -o ranked_fragments.json
```

