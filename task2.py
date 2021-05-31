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

"2. A dendrogram/heatmap representing the distance (global score) between ensembles."
def compute_global_score(med_mat_1, med_mat_2):
    ens_dRSMD = np.sqrt((1/N)*np.sum(np.power(med_mat_1 - med_mat_2,2))) #FORMULA 4 LAZAR
    return ens_dRSMD

def heatmap(pdb_ids, med_mats):
    heatmap = np.zeros((len(pdb_ids), len(pdb_ids)))
    for ens1 in range(len(pdb_ids)):
        for ens2 in range(len(pdb_ids)):
            #print("computing score between: ", ens1, ens2 )
            heatmap[ens1][ens2] = compute_global_score(med_mats[ens1], med_mats[ens2])

    fig, ax = plt.subplots()
    im = ax.imshow(heatmap)

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(pdb_ids)))
    ax.set_yticks(np.arange(len(pdb_ids)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(pdb_ids)
    ax.set_yticklabels(pdb_ids)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(pdb_ids)):
        for j in range(len(pdb_ids)):
            text = ax.text(j, i, "%.2f" %heatmap[i, j],
                           ha="center", va="center", color="w")

    ax.set_title("Global score between ensembles")
    fig.tight_layout()
    ensemble_string = "_".join(pdb_ids)
    plt.savefig("features/task2_heatmap_{}.png".format(ensemble_string))

def main():
    pdb_ids = []
    global M
    median_distance_mats = []
    for file in feature_files:
        pdb_ids.append(os.path.basename(file)[:-34]) # get pdb_id from file name
        with open(file) as f:
            data = json.load(f)
            M = len(data.keys())
            rg,mrasa, md, stdev_d = get_ensemble_features(data)
            median_distance_mats.append(md)
    heatmap(pdb_ids, median_distance_mats)



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
