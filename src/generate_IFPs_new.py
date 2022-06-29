import os
import argparse
import json
import pandas as pd
from joblib import Parallel, delayed
from oddt import toolkit
from tqdm import tqdm
from oddt.fingerprints_new import InteractionFingerprintModified


def get_bool(arg):
    if arg == 'true':
        return True
    if arg == 'false':
        return False


def get_ifp(pdb, mol, strict):

    lig = next(toolkit.readfile('mol', mol))
    try:
        prot = next(toolkit.readfile('pdb', pdb))
        prot.protein = True
    except:
        print('error', pdb)
    try:
        ifp = InteractionFingerprintModified(lig, prot, strict)
    except:
        print('error', pdb)

    return ifp


def process_files(files, input_dir, target, type='prot'):
    processed = []
    for f in files:
        spl = f.split('/capped/')
        new_f = os.path.join(input_dir, spl[-1])
        if target == 'nsp13' or target == 'PARP14A' and type == 'prot':
            new_f = new_f.replace('.holo_minimised_nolig.pdb', '_Repair.pdb')
        processed.append(new_f)

    return processed


def main():
    """
    Generate IFPs for merge results. Writes to json file.

    E.g. python src/generate_IFPs_new.py -f data/fragment_merging/Mpro.csv -s true -o data/fragment_merging/ifps
    -t Mpro -d /home/swills/Oxford/data/fragment_merging/run_2/output/capped -O Mpro
    """
    parser = argparse.ArgumentParser(
        epilog="python src/generate_IFPs_new.py -f data/fragment_merging/DPP11.csv -s true -o data/fragment_merging/ifps"
    )
    parser.add_argument('-t', '--target')
    parser.add_argument('-d', '--input_dir')
    parser.add_argument('-f', '--csv_file', help='csv file containing pdb files, mol files, smiles')
    parser.add_argument('-s', '--strict', help='true or false', default='false')
    parser.add_argument('-S', '--simsearch_dir', required=False)
    parser.add_argument('-o', '--output_dir', help='where to save output files')
    parser.add_argument('-O', '--output_fname', help='name of output json')
    args = parser.parse_args()

    csv = args.csv_file
    df = pd.read_csv(csv)

    dir = args.input_dir

    if args.simsearch_dir:
        pdbs = process_files(df['nolig_files'].tolist(), args.simsearch_dir, args.target)
    else:
        pdbs = process_files(df['nolig_files'].tolist(), dir, args.target)
    mols = process_files(df['mol_files'].tolist(), dir, args.target, 'mols')
    smiles = df['smiles'].tolist()

    strict = get_bool(args.strict)
    ifps = Parallel(n_jobs=os.cpu_count(),
                    backend="multiprocessing")(
        delayed(get_ifp)(pdb, mol, strict)
        for pdb, mol in tqdm(zip(pdbs, mols))
    )

    output_d = {}

    for smi, ifp in zip(smiles, ifps):
        if smi not in list(output_d.keys()):
            output_d[smi] = ifp
        else:
            n_old = len(output_d[smi])
            n_new = len(ifp)
            if n_new > n_old:
                output_d[smi] = ifp

    write_d = {args.target: output_d}
    json_file = os.path.join(args.output_dir,
                             f"ifps_{args.output_fname}.json")
    with open(json_file, "w") as f:
        json.dump(write_d, f)


if __name__ == "__main__":
    main()
