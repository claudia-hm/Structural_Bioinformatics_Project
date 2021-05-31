#!/usr/bin/env python

'''0. Utile functions'''
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


''' 1. Single conformation features'''
def get_radius_of_gyration(structure):
    '''
    Calculates the Radius of Gyration (Rg) of a protein in Angstroms.
    Does not use mass and assume heavy atoms have the same mass.
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
  hse = HSExposureCB(structure)
  rASA = np.zeros(N)
  i=0
  for chain in structure:
    for residue in chain:
        try:
            rASA[i]=hse[(chain.id, residue.id)][0]
            i=i+1
        except:
            print("")
  max=np.max(rASA)
  min=np.min(rASA)
  rASA= 1 - (rASA-min)/(max-min)
  return rASA.astype(np.float64)

def get_secondary_structure(structure):
    rama_ss_ranges = [(-180, -180, 80, 60, 'E', 'blue'),
                      (-180, 50, 80, 130, 'E', 'blue'),
                      (-100, -180, 100, 60, 'P', 'green'),
                      (-100, 50, 100, 130, 'P', 'green'),
                      (-180, -120, 180, 170, 'H', 'red'),
                      (0, -180, 180, 360, 'L', 'yellow')]

    # Calculate PSI and PHI
    ppb = PPBuilder()  # PolyPeptideBuilder
    ss = ["" for x in range(N)]
    for chain in structure:
        for pp in ppb.build_peptides(chain):
            phi_psi = pp.get_phi_psi_list()  # [(phi_residue_1, psi_residue_1), ...]
            for i, residue in enumerate(pp):
                    # print(model, chain, i, residue, phi_psi[i])
                    # Convert radians to degrees and remove first and last value that are None
                if phi_psi[i][0] is not None and phi_psi[i][1] is not None:
                    for x, y, w, h, ss_c, color in rama_ss_ranges:
                        if x <= phi_psi[i][0] < x + w and y <= phi_psi[i][1] < y + h:
                            ss[i]=ss_c
                            break
    return ss

def get_distance_matrix(structure, seq_sep=0):

    # Calculate the distance matrix
    distances = []
    for chain in structure:
        for residue1 in chain:
            if residue1.id[0] == " ":  # Exclude hetero/water residues
                row = []
                for residue2 in chain:
                    if residue2.id[0] == " ":  # Exclude hetero/water residues
                        if abs(residue1.id[1] - residue2.id[1]) >= seq_sep:
                            row.append(residue1["CA"] - residue2["CA"])
                        else:
                            row.append(None)
                distances.append(row)

    return np.array(distances).astype(np.float64)


def get_single_conformation_features(structure):
    rg = get_radius_of_gyration(structure)
    rasa = get_relative_asas(structure)
    ss = get_secondary_structure(structure)
    dm = get_distance_matrix(structure)

    return rg, rasa, ss, dm
'''2. Clustering functions '''
def mseDistanceMatrix(matrix1, matrix2):
    array1 = np.asarray(matrix1)
    array2 = np.asarray(matrix2)
    difference_array = np.subtract(array1, array2)
    squared_array = np.square(difference_array)
    mse = squared_array.mean()
    return mse


def mserASA(asa1, asa2):
    array1 = np.asarray(asa1)
    array2 = np.asarray(asa2)
    '''
    difference_array = np.subtract(array1, array2)
    squared_array = np.square(difference_array)
    mse = squared_array.mean()
    '''
    return np.abs(array1.mean() - array2.mean())


def ssNumberOfDismatch(ss1, ss2):
    max = 1
    countE = 0
    countP = 0
    countH = 0
    countL = 0
    countNone = 0
    np.count_nonzero(ss1 == ss2)
    for i in range(N):
        if ss1[i] == ss2[i]:
            if ss1[i] == "E":
                countE += 1
            elif ss1[i] == "P":
                countP += 1
            elif ss1[i] == "H":
                countH += 1
            elif ss1[i] == "L":
                countL += 1
            else:
                None
    return (N - (countE + countP + countH + countL)) / N


def getResidueRGdistanceMatrix():
    dgMat = np.zeros((M, M))
    for i in range(M):
        id = "{}_{}".format(pdb_id, i)
        for j in range(i, M):
            jd = "{}_{}".format(pdb_id, j)
            dgMat[i][j] = np.abs(np.subtract(features[id]['radius_of_giration'], features[jd]['radius_of_giration']))
    max = np.max(dgMat)
    min = np.min(dgMat)
    return np.divide(dgMat - min, max - min)


def getResidueASAdistanceMatrix():
    asaMat = np.zeros((M, M))
    for i in range(M):
        id = "{}_{}".format(pdb_id, i)
        for j in range(i, M):
            jd = "{}_{}".format(pdb_id, j)
            asaMat[i][j] = mserASA(features[id]['relative_asas'], features[jd]['relative_asas'])
    max = np.max(asaMat)
    min = np.min(asaMat)
    return np.divide(asaMat - min, max - min)


def getResidueSSdistanceMatrix():
    ssMat = np.zeros((M, M))
    for i in range(M):
        id = "{}_{}".format(pdb_id, i)
        for j in range(i, M):
            jd = "{}_{}".format(pdb_id, j)
            ssMat[i][j] = ssNumberOfDismatch(features[id]['secondary_structure'], features[jd]['secondary_structure'])
    max = np.max(ssMat)
    min = np.min(ssMat)
    return np.divide(ssMat - min, max - min)


def getResidueDMdistanceMatrix():
    distMat = np.zeros((M, M))
    for i in range(M):
        id = "{}_{}".format(pdb_id, i)
        for j in range(i, M):
            jd = "{}_{}".format(pdb_id, j)
            distMat[i][j] = mseDistanceMatrix(features[id]['distance_matrix'], features[jd]['distance_matrix'])
    max = np.max(distMat)
    min = np.min(distMat)
    return np.divide(distMat, max)

def get_distance_matrix_clustering():
    print("calcolo RG mat")
    dg = getResidueRGdistanceMatrix()

    print("calcolo ASA mat")
    da = getResidueASAdistanceMatrix()

    print("calcolo SS mat")
    ds = getResidueSSdistanceMatrix()

    print("calcolo DM mat")
    dm = getResidueDMdistanceMatrix()

    distanceMatrix = dg + da + ds + dm
    return distanceMatrix

def main():
    global pdb_id
    pdb_id = os.path.basename(pdb_path)[:-4]
    ensemble = PDBParser(QUIET=True).get_structure(pdb_id, pdb_path)
    print(pdb_id)
    global N
    global M
    N,M = get_ensemble_dimension(ensemble)
    rg = np.zeros(M)
    global features

    logging.info("Radius of gyration")
    for i in range(M):
        id = "{}_{}".format(pdb_id,i)
        features[id] = {}
        features[id]['radius_of_giration'] = get_radius_of_gyration(ensemble[i])
        features[id]['relative_asas'] = get_relative_asas(ensemble[i]).tolist()
        features[id]['secondary_structure'] = list(get_secondary_structure(ensemble[i]))
        features[id]['distance_matrix'] =get_distance_matrix(ensemble[i]).tolist()



    with open("{}_single_conformation_features.json".format(pdb_id),'w') as outfile:
        json.dump(features, outfile)

    import matplotlib.pyplot as plt

    clustering_dm = get_distance_matrix_clustering()
    plt.imshow(clustering_dm, cmap='hot', interpolation='nearest')
    plt.show()

if __name__ == "__main__":
    import argparse
    import logging
    import json
    import numpy as np
    from Bio.PDB import is_aa, PPBuilder, HSExposureCB
    from Bio.PDB.PDBParser import PDBParser
    import math
    import argparse
    import os

    N,M = 0,0
    features = {}
    pdb_id = ""
    logging.basicConfig(filename="task1.log", encoding='utf-8', level=logging.INFO,
                        format="%(asctime)s %(levelname)s: %(message)s", filemode='w')

    parser = argparse.ArgumentParser(description='Code for task 1 Structural Bioinformatics Project')
    parser.add_argument('file', metavar='F', nargs=1,
                        help='PDB file path')

    args = parser.parse_args()
    pdb_path = args.file[0]
    logging.info("Program start")
    main()