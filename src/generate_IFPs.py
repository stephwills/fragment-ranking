import os
from oddt.fingerprints_new import InteractionFingerprint, InteractionFingerprintAtomic, tanimoto
import oddt
import sys 
from pymol import cmd
import statistics
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem import DataStructs, SaltRemover
import json
import argparse

def separate_files(filepath):

    # takes in a protein-ligand complex and separates files into protein and ligand pdb files

    cmd.reinitialize()
    cmd.load(filepath, 'complex')
    cmd.select('lig', 'resn LIG')
    cmd.save('data/tmp/lig.pdb', 'lig')
    cmd.extract('hets', 'complex and HETATM')
    cmd.save('data/tmp/prot.pdb', 'complex')


def get_IFP(IFP_type):

    # uses the previously split protein and ligand files and caculates binary protein-ligand interaction fingerprint

    lig = next(oddt.toolkit.readfile('pdb', 'data/tmp/lig.pdb'))
    prot = next(oddt.toolkit.readfile('pdb', 'data/tmp/prot.pdb'))

    prot.protein = True

    if IFP_type == 'residue':
        IFP = InteractionFingerprint(lig, prot, strict=False)
    elif IFP_type == 'atomic':
        IFP = InteractionFingerprintAtomic(lig, prot, strict=False)

    return IFP


def get_library_smiles(library_sdf):

    # using the sdf file of the fragment library, the molecules are cleaned of salts and converted into SMILES strings

    remover = SaltRemover.SaltRemover()

    raw_mols = list(Chem.SDMolSupplier(library_sdf))
    clean_smiles = [remover.StripMol(mol) for mol in raw_mols]
    library_smiles = [Chem.MolToSmiles(i) for i in clean_smiles]

    return library_smiles


def get_IFP_vectors(data_path, target, IFP_type):

    # all IFPs for a particular target are calculated

    ismiles = []
    ifrags = []
    ivecs = []

    for ligand in os.listdir(f'{data_path}/{target}'):


        try:
            if Chem.MolToSmiles(Chem.MolFromSmiles(xtal_smiles[ligand])) in library_smiles:

                separate_files(f'{data_path}/{target}/{ligand}')

                IFP = get_IFP(IFP_type)
            
                if list(IFP).count(0) < len(IFP):
                    ismiles.append(Chem.MolToSmiles(Chem.MolFromSmiles(xtal_smiles[ligand])))
                    ifrags.append(ligand)
                    ivecs.append(IFP)
                    print(ligand, [i for i in range(len(IFP)) if IFP[i] > 0])
                else:
                    print(ligand, 'no interactions detected')

        
        except:
            pass
    
    return ismiles, ifrags, ivecs


def get_uniform_IFPs(ismiles, ifrags, ivecs):

    # IFPs are analysed to determine most common length (due to some models missing residues)
    # only IFPs of identical length are returned

    smiles = []
    frags = []
    vecs = []

    lengths = [len(i) for i in ivecs]
    length = statistics.mode(lengths)

    wrong = 0
    for i in range(len(ivecs)):
        if len(ivecs[i]) == length:
            vecs.append(ivecs[i])
            frags.append(ifrags[i])
            smiles.append(ismiles[i])
        else:
            wrong += 1

    return vecs, frags, smiles, wrong


def get_smiles_bits(vecs, smiles):

    # take IFP vectors and assign 'on' bits to the smiles strings responsible for the interaction

    smiles_bits = {}

    for i in range(len(smiles)):
        if smiles[i] not in smiles_bits:
            smiles_bits[smiles[i]] = []

        for f in range(len(smiles)):
            mol1, mol2 = Chem.MolFromSmiles(smiles[i]), Chem.MolFromSmiles(smiles[f])
            mol_sim = DataStructs.DiceSimilarity(AllChem.GetMorganFingerprint(mol1,2), AllChem.GetMorganFingerprint(mol2, 2))
            if mol_sim == 1:
                smiles_bits[smiles[i]] += [b for b in range(len(vecs[f])) if vecs[f][b] != 0]
        
    for smiles in smiles_bits:
        smiles_bits[smiles] = list(set(smiles_bits[smiles]))
    
    return smiles_bits


if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='atomic or residue interactions?')


    parser.add_argument('-IFP', '--ifp')
    parser.add_argument('-sdf', '--SDF')
    parser.add_argument('-exps', '--EXPS')
    parser.add_argument('-pdbs', '--PDBS')

    args = vars(parser.parse_args())

    IFP_type = args['ifp']
    library_sdf = args['SDF']
    experiments_json = args['EXPS']
    data_path = args['PDBS']

    assert IFP_type == 'residue' or IFP_type == 'atomic'

    target_data = {}
    library_smiles = get_library_smiles(library_sdf)
    print('Total unique compounds in library:', len(set(library_smiles)))

    xtal_smiles = json.load(open(experiments_json, 'r'))


    for target in os.listdir(data_path):

        try:
            print(target)

            ifrags, ivecs, ismiles = get_IFP_vectors(data_path, target, IFP_type)
            vecs, frags, smiles, wrong = get_uniform_IFPs(ifrags, ivecs, ismiles)

            print('Useable IFPs:', len(vecs))
            print('Structures with missing atoms:', wrong)

            smiles_bits = get_smiles_bits(vecs, smiles)
            target_data[target] = smiles_bits

            if IFP_type == 'residue':
                json.dump(target_data, open('data/smiles_bits_residue.json', 'w'))
            elif IFP_type == 'atomic':
                json.dump(target_data, open('data/smiles_bits_atomic.json', 'w'))
        
            print(f'{target} complete')
        except:
            print(target, 'error')
            print(sys.exc_info()[1])

