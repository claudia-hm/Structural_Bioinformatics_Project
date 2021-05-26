#!/usr/bin/env python
def get_radius_of_gyration(structure):
    '''
    Calculates the Radius of Gyration (Rg) of a protein in Angstroms.
    Does not use mass and assume heavy atoms have the same mass.

    https://en.wikipedia.org/wiki/Radius_of_gyration  (formula considering mass)
    https://link.springer.com/article/10.1134/S0026893308040195  (formula without mass)
    '''

    # Heavy atoms coordinates
    coord = list()
    for atom in structure.get_atoms():
        if atom.get_name()[0] in ['C', 'O', 'N', 'S']:
            coord.append(atom.get_coord())
    coord = np.array(coord)  # N X 3

    barycenter = np.sum(coord, axis=0) / coord.shape[0]  # center of mass is more correct

    # Calculate distance of each atom from the barycenter
    dist = coord - barycenter
    dist = dist * dist
    dist = np.sqrt(np.sum(dist, axis=1))
    #print(dist)

    return round(math.sqrt(np.sum(dist * dist) / len(coord)), 3)

def get_relative_asas(structure):
    rASA = np.zeros(N)
    dssp = DSSP(structure, pdb_path, dssp="binx/dssp-2.3.0/mkdssp")  # WARNING Check the path of mkdssp
    i=0
    for ss in dssp:
        rASA[i]=ss[3]
        i = i+1
    return rASA

def get_relative_asas2(structure):
  hse = HSExposureCB(structure)
  rASA = np.zeros(N)
  i=0
  for chain in structure:
    for residue in chain:
        try:
            rASA[i]=hse[(chain.id, residue.id)][0]
            i=i+1
        except:
            print(residue)
  max=np.max(rASA)
  min=np.min(rASA)
  rASA= 1 - (rASA-min)/(max-min)
  return rASA

def get_secondary_structure(structure):
    rama_ss_ranges = [(-180, -180, 80, 60, 'E', 'blue'),
                      (-180, 50, 80, 130, 'E', 'blue'),
                      (-100, -180, 100, 60, 'P', 'green'),
                      (-100, 50, 100, 130, 'P', 'green'),
                      (-180, -120, 180, 170, 'H', 'red'),
                      (0, -180, 180, 360, 'L', 'yellow')]

    num_struct = N

    # Calculate PSI and PHI
    ppb = PPBuilder()  # PolyPeptideBuilder
    ss = []

    for chain in structure:
        for pp in ppb.build_peptides(chain):
            phi_psi = pp.get_phi_psi_list()  # [(phi_residue_1, psi_residue_1), ...]
            for i, residue in enumerate(pp):
                    # print(model, chain, i, residue, phi_psi[i])
                    # Convert radians to degrees and remove first and last value that are None
                if phi_psi[i][0] is not None and phi_psi[i][1] is not None:
                    for x, y, w, h, ss_c, color in rama_ss_ranges:
                        if x <= phi_psi[i][0] < x + w and y <= phi_psi[i][1] < y + h:
                            ss.append(ss_c)
                            break

    return ss

def get_distance_matrix(structure, seq_sep=6):

    # Calculate the distance matrix
    distances = []
    for residue1 in structure:
        if residue1.id[0] == " ":  # Exclude hetero/water residues
            row = []
            for residue2 in structure:
                if residue2.id[0] == " ":  # Exclude hetero/water residues
                    if abs(residue1.id[1] - residue2.id[1]) >= seq_sep:
                        row.append(residue1["CA"] - residue2["CA"])
                    else:
                        row.append(None)
            distances.append(row)

    return distances


def get_single_conformation_features(structure):
    rg = get_radius_of_gyration(structure)
    rasa = get_relative_asas(structure)
    ss = get_secondary_structure(structure)
    dm = get_distance_matrix(structure)

    return rg, rasa, ss, dm

def get_ensemble_dimension(ensemble):
    structure = ensemble[0]
    N,M = 0,len(ensemble)
    for chain in structure:
        for residue in chain:
            if is_aa(residue):  # Filter hetero groups (returns only amino acids)                                      IUPACData.protein_letters_3to1.get(residue.get_resname().capitalize())))
                N += 1
            else:
                print(residue.id)
    return N,M

def main():
    pdb_id = os.path.basename(pdb_path)[:-4]
    ensemble = PDBParser(QUIET=True).get_structure(pdb_id, pdb_path)
    logging.basicConfig(filename="task1.log", encoding ='utf-8', level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s", filemode='w')
    # logging.debug(pdb_id)
    # logging.warning("aaaaah")
    logging.info("Program start")
    # logging.error("bbbbbb")
    global N
    global M
    N,M = get_ensemble_dimension(ensemble)
    rg = np.zeros(M)
    features = {}
    logging.info("Radius of gyration")
    for i in range(M):
        id = "{}_{}".format(pdb_id,i)
        features[id] = {}
        features[id]['radius_of_giration'] = get_radius_of_gyration(ensemble[i])
        features[id]['relative_asas'] = list(get_relative_asas(ensemble[i]))
        features[id]['secondary_structure'] = list(get_secondary_structure(ensemble[i]))
        features[id]['distance_matrix'] = get_distance_matrix(ensemble[i])



    with open("{}_single_conformation_features.txt".format(pdb_id),'w') as outfile:
        json.dump(features, outfile)




if __name__ == "__main__":
    import argparse
    import logging
    import json
    import numpy as np
    from Bio.PDB import PDBList, Superimposer, is_aa, PPBuilder, HSExposureCB
    from Bio.PDB.DSSP import DSSP
    from Bio.PDB.PDBParser import PDBParser
    import requests
    import math
    import matplotlib.pyplot as plt
    import argparse
    import os

    N,M = 0,0

    parser = argparse.ArgumentParser(description='Code for task 1 Structural Bioinformatics Project')
    parser.add_argument('file', metavar='F', nargs=1,
                        help='PDB file path')

    args = parser.parse_args()
    pdb_path = args.file[0]
    main()