# fragment-ranking

A method for ranking fragments by how much novel information they give about protein targets in fragment screens. When using the results of fragment screens on many diverse targets, this method has been shown to select a set of functionally diverse fragments that can get information more efficiently from new targets. 

A custom version of [ODDT](https://github.com/oddt/oddt)[1] is included, which contains a module that generates an atomic-level interaction fingerprint.

## Instructions for use

Conda or miniconda is required. Downloards and instructions can be found at: https://docs.conda.io/projects/conda/en/latest/. Once this requirement has been fulfilled, you can install fragment-ranking. 

Installation:

1. clone repository
```bash
git clone https://github.com/annacarbery/fragment-ranking.git
cd fragment-ranking
```
2. create new conda environment and install pymol and rdkit
```bash
conda create -n fragment-ranking -y
conda activate fragment-ranking
conda install -c schrodinger pymol -y
conda install -c conda-forge rdkit -y 
```
3. Install custom ODDT version (provided within fragment-ranking repository)
```bash
cd oddt
python setup.py install
cd ..
```

Required input files:
1. Fragment library, in sdf format.
2. Structures of targets with fragments bound, in pdb format. Structures should not have any atoms or residues missing.
3. JSON dictionary of pdb file with associated fragment (in SMILES string), for example : {"mArh-x1018.pdb": "Clc1ccc2nnnn2n1",...}

Test data is provided.

Run:
1. Generate interaction fingerprints from structures. Either residue-level or atomic-level interaction fingerprints can be used. 
```bash
python src/generate_IFPs.py -IFP [atomic/residue] -sdf [path_to_fragment_library.sdf] -exps [path_to_experiments.json] -pdbs [path_to_pdb_files]
```
  for example:
```bash
python src/generate_IFPs.py -IFP atomic -sdf data/library.sdf -exps data/experiments_mArh.json -pdbs data/structures/
```
2. Rank fragments
```bash
python src/rank_fragments.py -sdf [path_to_fragment_library.sdf] -bits data/smiles_bits_[atomic/residue].json -o [output_file]
```
for example:
```bash
python src/rank_fragments.py -sdf data/library.sdf -bits data/smiles_bits_atomic.json -o ranked_fragments.json
```

# References

1. Wójcikowski, M., Zielenkiewicz, P., & Siedlecki, P. (2015). Open Drug Discovery Toolkit (ODDT): a new open-source player in the drug discovery field. Journal of Cheminformatics, 7(1), 26. doi:10.1186/s13321-015-0078-2
