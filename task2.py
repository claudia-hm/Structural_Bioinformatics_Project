#!/usr/bin/env python

'''0. Utile functions'''




'''1. Ensembles features (multiple conformations)'''
def get_radius_of_gyration(data):
    global M
    radius_gyration_vec = np.zeros(M)
    for i, key in enumerate(data):
        radius_gyration_vec[i] = data[key]["radius_of_giration"]
    return radius_gyration_vec

def get_median_asa(data):
    global N
    N = len(data[list(data.keys())[0]]['relative_asas'])

    med_rasa = np.zeros(N)
    for i in range(N):  # for
        rasas = []
        for key in data:
            rasas.append(data[key]["relative_asas"][i])
        med_rasa[i]= np.median(rasas)
    return med_rasa

def get_median_distance(data):
    global N
    N = len(data[list(data.keys())[0]]['distance_matrix'])

    med_distance = np.zeros((N, N))
    for i in range(N):  # for
        for j in range(N):  # for every conformation
            distances = []
            for key in data:
                distances.append(data[key]["distance_matrix"][i][j])
            med_distance[i, j] = np.median(distances)
    return med_distance



def get_stdev_distance(data):
    global N
    N = len(data[list(data.keys())[0]]['distance_matrix'])

    stdev_distance = np.zeros((N, N))
    for i in range(N):  # for
        for j in range(N):  # for every conformation
            distances = []
            for key in data:
                distances.append(data[key]["distance_matrix"][i][j])
            stdev_distance[i, j] = np.std(distances)
    return stdev_distance

def get_ensemble_features(data):

    rg = get_radius_of_gyration(data) # 1. Radius of gyration for each conformation in the ensemble.
    mrasa = get_median_asa(data) #3. Median solvent accessibility for each position across ensemble conformations.
    md = get_median_distance(data)  # 5. Median distance of each pair of equivalent positions across ensemble conformations.
    stdev_d = get_stdev_distance(data)  # 6. Standard deviation of the distance of each pair of equivalent positions across ensemble conformations.

    return rg,mrasa, md, stdev_d

def main():
    pdb_ids = []
    global M
    for file in feature_files:
        pdb_ids.append(os.path.basename(file)[:-34]) # get pdb_id from file name
        with open(file) as f:
            data = json.load(f)
            M = len(data.keys())
            get_ensemble_features(data)





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
    logging.basicConfig(filename="task2.log", encoding ='utf-8', level=logging.INFO, format="%(asctime)s %(levelname)s: %(message)s", filemode='w')

    parser = argparse.ArgumentParser(description='Code for task 2 Structural Bioinformatics Project')
    parser.add_argument('feature_files', metavar='F', nargs='+',
                        help='Feature files')
    args = parser.parse_args()
    if len(args.feature_files) > 1:
        feature_files = args.feature_files
        print(feature_files)
        logging.info("Program start")
        main()
    else:
        logging.error("Not enough number of feature files (minimum required is 2)")
