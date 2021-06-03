#!/usr/bin/env python

'''0. Utile functions'''



'''1. Ensembles features (multiple conformations)'''
def get_radius_of_gyration(data):
    global M
    radius_gyration_vec = np.zeros(M)
    for i, key in enumerate(data):
        radius_gyration_vec[i] = data[key]["radius_of_giration"]
    return radius_gyration_vec

def get_ss_entropy(data):

    # input matrix (MxN): M ensemble conformation, N residue
    # list=[E,P,H,L,nn,entropy]
    global N
    N = len(data[list(data.keys())[0]]["secondary_structure"])

    ss_mat = []
    for key in data:
        ss_mat.append(data[key]["secondary_structure"])

    num_row, num_col = M, N

    output = np.zeros((num_col, 6))

    for int_mod in range(0, num_row):
        for int_residue in range(0, num_col):

            # ss_residue=data[int_mod,int_residue]
            ss_residue = ss_mat[int_mod][int_residue]
            # print(type(ss_residue))

            if ss_residue == 'E':
                output[int_residue][0] += 1
            elif ss_residue == 'P':
                output[int_residue][1] += 1
            elif ss_residue == 'H':
                output[int_residue][2] += 1
            elif ss_residue == 'L':
                output[int_residue][3] += 1
            else:
                output[int_residue][4] += 1

    for int_residue in range(0, num_col):

        output[int_residue] = output[int_residue] / num_row

        log_frq = np.array([1, 1, 1, 1, 1], dtype=np.float)

        for i in range(5):
            if output[int_residue, i] != 0:
                log_frq[i] = np.log(output[int_residue, i])
        output[int_residue, 5] = np.sum(-output[int_residue, 0:5] * log_frq)

    return output[:, 5]

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

def get_median_rsmd(data, pdb_id):
    ensemble = PDBParser(QUIET=True).get_structure(pdb_id, "data/{}.pdb".format(pdb_id))
    super_imposer = Superimposer()
    window_size = 5

    _RMSD_along_all_conf = []  # list of RMSD_n values acrossall possibile pairs of conformation

    for j, model_j in enumerate(ensemble):
        if j > 0:
            model_rmsd = []  # RMSD, no_fragment X fragment_size
            alt_model = [atom for atom in model_j.get_atoms() if atom.get_name() == "CA"]  # coords of the model

            # Iterate fragments
            for start in range(len(ref_model) - window_size):
                end = start + window_size
                ref_fragment = ref_model[start:end]
                alt_fragment = alt_model[start:end]

                # Calculate rotation/translation matrices
                super_imposer.set_atoms(ref_fragment, alt_fragment)
                # print(super_imposer.rms, super_imposer.rotran)

                # Rotate-translate coordinates
                alt_fragment_coord = np.array([atom.get_coord() for atom in alt_fragment])
                alt_fragment_coord = np.dot(super_imposer.rotran[0].T, alt_fragment_coord.T).T
                alt_fragment_coord = alt_fragment_coord + super_imposer.rotran[1]

                # Calculate RMSD
                # https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
                ref_fragment_coord = np.array([atom.get_coord() for atom in ref_fragment])
                dist = ref_fragment_coord - alt_fragment_coord
                rmsd_fragment = np.sqrt(
                    np.sum(dist * dist) / window_size)  # Total RMSD of the fragment. Identical to super_imposer.rms
                # rmsd_res = np.sqrt(np.sum(dist * dist, axis=1))  # RMSD for each residue of the fragment
                model_rmsd.append(rmsd_fragment)
            # print ("modeli-{}, modelj-{},".format(model_i,model_j))
            _RMSD_along_all_conf.append(model_rmsd)
        else:
            ref_model = [atom for atom in model_j.get_atoms() if atom.get_name() == "CA"]  # CA of the first model
    np_RMSD_along_all_conf_easy = np.array(_RMSD_along_all_conf)

    num_residue = np_RMSD_along_all_conf_easy.shape[1]

    mean_RMSD_along_all_conf_easy = np.mean(np_RMSD_along_all_conf_easy, axis=0)
    median_RMSD_along_all_conf_easy = np.median(np_RMSD_along_all_conf_easy, axis=0)
    stDev_RMSD_along_all_conf_easy = np.std(np_RMSD_along_all_conf_easy, axis=0)
    ci_easy = 1.96 * stDev_RMSD_along_all_conf_easy / mean_RMSD_along_all_conf_easy
    p90_RMSD_along_all_conf_easy = np.percentile(np_RMSD_along_all_conf_easy, [5, 95], axis=0)
    p95_RMSD_along_all_conf_easy = np.percentile(np_RMSD_along_all_conf_easy, [2.5, 97.5], axis=0)
    return median_RMSD_along_all_conf_easy

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




def get_ensemble_features(data, pdb_id):

    rg = get_radius_of_gyration(data) # 1. Radius of gyration for each conformation in the ensemble.
    ss = get_ss_entropy(data)# 2. Secondary structure entropy for each position across ensemble conformations.
    mrasa = get_median_asa(data) #3. Median solvent accessibility for each position across ensemble conformations.
    m_rsmd = get_median_rsmd(data, pdb_id)# 4. Median RMSD for each position across ensemble conformations 1
    md = get_median_distance(data)  # 5. Median distance of each pair of equivalent positions across ensemble conformations.
    stdev_d = get_stdev_distance(data)  # 6. Standard deviation of the distance of each pair of equivalent positions across ensemble conformations.

    return rg, ss, mrasa, m_rsmd, md, stdev_d

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
        pdb_id = os.path.basename(file)[:-34]
        pdb_ids.append(pdb_id) # get pdb_id from file name
        with open(file) as f:
            data = json.load(f)
            M = len(data.keys())
            rg, ss, mrasa, m_rsmd, md, stdev_d = get_ensemble_features(data, pdb_id)
            median_distance_mats.append(md)
        features = {}

        features['radius_of_giration'] = rg.tolist()
        features['ss_entropy'] = ss.tolist()
        features['median_solvent_accessibility'] = mrasa.tolist()
        features['median_rsmd'] = m_rsmd.tolist()
        features['median_distance'] = md.tolist()
        features['stdev_distance'] = stdev_d.tolist()

        with open("features/{}_ensemble_features.json".format(pdb_id),'w') as outfile:
            json.dump(features, outfile)




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
