import json
from tqdm import tqdm
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs, SaltRemover
import os
from generate_IFPs import get_library_smiles
import argparse


def get_next_best_compound(smiles_bits, current_info, current_compounds, all_compounds):
    # for each library size, selects new compound that gives most additional information

    new_compounds = [i for i in all_compounds if i not in current_compounds]  # current compounds initially empty

    best_new_info = 0
    best_new_compound = None
    best_new_bits = {}

    for compound in new_compounds:  # in compounds not ranked yet
        new_info = 0

        for target in smiles_bits:
            if compound in smiles_bits[target]:
                # get number interactions made by compound if interaction hasn't already been made by another compound
                new_info += len([i for i in smiles_bits[target][compound] if i not in current_info[target]])

        if new_info > best_new_info:
            # if compound makes largest number of new interactions so far, this becomes the best one
            best_new_info = new_info
            best_new_compound = compound
            best_new_bits = {}  # best new bits recorded by target {target: [bits for that target...]}
            for target in smiles_bits:
                if compound in smiles_bits[target]:
                    best_new_bits[target] = [i for i in smiles_bits[target][compound] if i not in current_info[target]]

    for target in best_new_bits:
        current_info[target] += best_new_bits[target]  # add the interactions made to the info recorded so far

    return current_info, best_new_compound


def rank(smiles_bits, all_smiles, output):
    # ranks fragments provided by how much novel information they give across all targets

    fraction = [0]
    ranked_frags = []

    current_info = {}
    for t in smiles_bits:
        current_info[t] = []  # {target: [] ...}

    for s in tqdm(range(len(all_smiles))):
        current_info, best_new_compound = get_next_best_compound(smiles_bits, current_info, ranked_frags, all_smiles)
        ranked_frags.append(best_new_compound)

    json.dump(ranked_frags, open(output, 'w'))

    return ranked_frags


if __name__ == '__main__':
    parser = argparse.ArgumentParser()

    parser.add_argument('-i', '--ifp_file')
    parser.add_argument('-o', '--output')

    args = vars(parser.parse_args())

    smiles_bits = args['ifp_file']
    smiles_bits = json.load(open(smiles_bits, 'r'))

    library_smiles = []
    for target in smiles_bits:
        for smi in smiles_bits[target]:
            library_smiles.append(smi)

    output = args['output']
    ranked_fragments = rank(smiles_bits, library_smiles, output)

    print('fragments ranked')
