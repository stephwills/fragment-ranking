import json
import matplotlib.pyplot as plt
import numpy as np 
import random
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs, SaltRemover
import os
from generate_IFPs import get_library_smiles\
import argparse

def get_next_best_compound(smiles_bits, current_info, current_compounds, all_compounds):

    # for each library size, selects new compound that gives most additional information
    
    new_compounds = [i for i in all_compounds if i not in current_compounds]

    best_new_info = 0
    best_new_compound = None
    best_new_bits = {}

    for compound in new_compounds:
        new_info = 0

        for target in smiles_bits:
            if compound in smiles_bits[target]:
                new_info += len([i for i in smiles_bits[target][compound] if i not in current_info[target]])

        if new_info > best_new_info:
            best_new_info = new_info
            best_new_compound = compound
            best_new_bits = {}
            for target in smiles_bits:
                if compound in smiles_bits[target]:
                    best_new_bits[target] = [i for i in smiles_bits[target][compound] if i not in current_info[target]]

        
    for target in best_new_bits:
        current_info[target] += best_new_bits[target]

    return current_info, best_new_compound        



def rank(smiles_bits, all_smiles, output):

    # ranks fragments provided by how much novel information they give across all targets

    fraction = [0]
    comps = []

    current_info = {}
    for t in smiles_bits:
        current_info[t] = []

    for s in range(len(all_smiles)):

        current_info, best_new_compound = get_next_best_compound(smiles_bits, current_info, ranked_frags, all_smiles)
        ranked_frags.append(best_new_compound)

    json.dump(ranked_frags, open(output, 'w'))

    return ranked_frags


if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    parser.add_argument('-sdf', '--SDF')
    parser.add_argument('-bits', '--BITS')
    parser.add_argument('-o', '--output')

    args = vars(parser.parse_args())

    library_sdf = args['SDF']
    smiles_bits = args['BITS']
    output = args['output']

    library_smiles = get_library_smiles(library_sdf)
    ranked_fragments = rank(smiles_bits, library_smiles, output)

    print('fragments ranked')
