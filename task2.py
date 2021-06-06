#!/usr/bin/env python


'''1. Ensembles features (multiple conformations)'''
def get_radius_of_gyration(data):
    logging.info("Computing radius of gyration")
    radius_gyration_vec = np.zeros(M)
    for i, key in enumerate(data):
        radius_gyration_vec[i] = data[key]["radius_of_giration"]
    return radius_gyration_vec

def get_ss_entropy(data):
    logging.info("Computing Secondary Structure entropy")
    global N
    N = len(data[list(data.keys())[0]]["secondary_structure"])
    ss_mat = []
    for key in data:
        ss_mat.append(data[key]["secondary_structure"])
    num_row, num_col = M, N
    output = np.zeros((num_col, 6))
    for int_mod in range(0, num_row):
        for int_residue in range(0, num_col):
            ss_residue = ss_mat[int_mod][int_residue]
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
    logging.info("Computing median solvent accessibility ")
    global N
    N = len(data[list(data.keys())[0]]['relative_asas'])
    med_rasa = np.zeros(N)
    for i in range(N):  # for
        rasas = []
        for key in data:
            rasas.append(data[key]["relative_asas"][i])
        med_rasa[i]= np.median(rasas)
    return med_rasa

def get_median_rsmd(pdb_id):
    logging.info("Computing median RSMD ")
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

    mean_RMSD_along_all_conf_easy = np.mean(np_RMSD_along_all_conf_easy, axis=0)
    median_RMSD_along_all_conf_easy = np.median(np_RMSD_along_all_conf_easy, axis=0)
    stDev_RMSD_along_all_conf_easy = np.std(np_RMSD_along_all_conf_easy, axis=0)
    return median_RMSD_along_all_conf_easy

def get_median_distance(data):
    logging.info("Computing median distance")
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
    logging.info("Computing standard deviation distance")
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
    m_rsmd = get_median_rsmd(pdb_id)# 4. Median RMSD for each position across ensemble conformations 1
    md = get_median_distance(data)  # 5. Median distance of each pair of equivalent positions across ensemble conformations.
    stdev_d = get_stdev_distance(data)  # 6. Standard deviation of the distance of each pair of equivalent positions across ensemble conformations.

    return rg, ss, mrasa, m_rsmd, md, stdev_d

"2. A dendrogram/heatmap representing the distance (global score) between ensembles."
def heatmap(pdb_ids, global_score):
    logging.info("Creating global score heatmap")
    fig, ax = plt.subplots()
    im = ax.imshow(global_score)

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
            text = ax.text(j, i, "%.2f" %global_score[i, j],
                           ha="center", va="center", color="w")

    ax.set_title("Global score between ensembles")
    fig.tight_layout()
    ensemble_string = "_".join(pdb_ids)
    plt.savefig("output/task2_heatmap_{}.png".format(ensemble_string))


def compute_global_score(med_mat_1, med_mat_2):
    ens_dRSMD = np.sqrt((1/N)*np.sum(np.power(np.array(med_mat_1) - np.array(med_mat_2),2))) #FORMULA 4 LAZAR
    return ens_dRSMD

def get_global_score(pdb_ids, ensemble_features):
    global_score = np.zeros((len(pdb_ids), len(pdb_ids)))
    for ens1 in range(len(pdb_ids)):
        for ens2 in range(len(pdb_ids)):
            logging.info("Computing global score between ensemble")
            global_score[ens1][ens2] = compute_global_score(ensemble_features[ens1]["median_distance"], ensemble_features[ens2]["median_distance"])
    return global_score

def get_RMSD_local_score(pdb_ids, ensembles_features):
    rmsd_local_score = {}
    for ens1 in range(len(pdb_ids)):
        for ens2 in range(ens1+1,len(pdb_ids)):
            rmsd1 = ensembles_features[ens1]['median_rsmd']
            rmsd2 = ensembles_features[ens2]['median_rsmd']
            rmsd_local_score[(pdb_ids[ens1], pdb_ids[ens2])] = np.abs(np.subtract(rmsd1,rmsd2))
    return rmsd_local_score

def get_rASA_local_score(pdb_ids, ensembles_features):
    rmsd_local_score = {}
    for ens1 in range(len(pdb_ids)):
        for ens2 in range(ens1+1,len(pdb_ids)):
            rmsd1 = ensembles_features[ens1]['median_solvent_accessibility']
            rmsd2 = ensembles_features[ens2]['median_solvent_accessibility']
            rmsd_local_score[(pdb_ids[ens1], pdb_ids[ens2])] = np.abs(np.subtract(rmsd1,rmsd2))
    return rmsd_local_score

def get_ss_entropy_local_score(pdb_ids, ensembles_features):
    rmsd_local_score = {}
    for ens1 in range(len(pdb_ids)):
        for ens2 in range(ens1+1,len(pdb_ids)):
            rmsd1 = ensembles_features[ens1]['ss_entropy']
            rmsd2 = ensembles_features[ens2]['ss_entropy']
            rmsd_local_score[(pdb_ids[ens1], pdb_ids[ens2])] = np.abs(np.subtract(rmsd1,rmsd2))
    return rmsd_local_score


def plot_local_score(pdb_ids, ensembles_features, ens1, ens2):
    logging.info("Computing RSMD local score")
    local_score_dic = get_RMSD_local_score(pdb_ids, ensembles_features)
    logging.info("Computing rASA local score")
    local_score_dic2 = get_rASA_local_score(pdb_ids, ensembles_features)
    logging.info("Computing SS entropy local score")
    local_score_dic3 = get_ss_entropy_local_score(pdb_ids, ensembles_features)
    local_score1 = local_score_dic[(pdb_ids[ens1], pdb_ids[ens2])]
    local_score2 = local_score_dic2[(pdb_ids[ens1], pdb_ids[ens2])]
    local_score3 = local_score_dic3[(pdb_ids[ens1], pdb_ids[ens2])]

    logging.info("Creating local score plot")
    median_rsmd1 = np.array(ensembles_features[ens1]['median_rsmd'])
    median_rsmd2 = np.array(ensembles_features[ens2]['median_rsmd'])

    rasa1 = np.array(ensembles_features[ens1]['median_solvent_accessibility'])
    rasa2 = np.array(ensembles_features[ens2]['median_solvent_accessibility'])

    entropy1 = np.array(ensembles_features[ens1]['ss_entropy'])
    entropy2 = np.array(ensembles_features[ens2]['ss_entropy'])

    fig, (ax1, ax2, ax3) = plt.subplots(3, 2, sharex=True)
    fig.suptitle("Local scores {}_{}".format(pdb_ids[ens1], pdb_ids[ens2]))
    fig.set_figheight(10)
    fig.set_figwidth(20)
    ax1[0].plot(range(N-5), median_rsmd2, color='#990000', linewidth=1, label=pdb_ids[ens1])
    ax1[0].plot(range(N-5), median_rsmd1, color='#000099', linewidth=1, label=pdb_ids[ens2])
    ax1[0].fill_between(range(N-5), median_rsmd1, median_rsmd2, where=median_rsmd2 <= median_rsmd1, facecolor='#dddddd',
                        interpolate=True)
    ax1[0].fill_between(range(N-5), median_rsmd1, median_rsmd2, where=median_rsmd2 >= median_rsmd1, facecolor='#dddddd',
                        interpolate=True)
    ax1[0].set_ylabel('RSMD local score')
    ax1[0].legend()

    ax1[1].plot(range(N-5), local_score1, color='black', linewidth=1)

    ax2[0].plot(range(N), rasa2, color='#990000', linewidth=1, label=pdb_ids[ens1])
    ax2[0].plot(range(N), rasa1, color='#000099', linewidth=1, label=pdb_ids[ens2])
    ax2[0].fill_between(range(N), rasa1, rasa2, where=rasa2 <= rasa1, facecolor='#dddddd', interpolate=True)
    ax2[0].fill_between(range(N), rasa1, rasa2, where=rasa2 >= rasa1, facecolor='#dddddd', interpolate=True)
    ax2[0].set_ylabel('Median rasa local score')
    ax2[0].legend()

    ax2[1].plot(range(N), local_score2, color='black', linewidth=1)

    ax3[0].plot(range(N), entropy2, color='#990000', linewidth=1, label=pdb_ids[ens1])
    ax3[0].plot(range(N), entropy1, color='#000099', linewidth=1, label=pdb_ids[ens2])
    ax3[0].fill_between(range(N), entropy1, entropy2, where=entropy2 <= entropy1, facecolor='#dddddd',
                        interpolate=True)
    ax3[0].fill_between(range(N), entropy1, entropy2, where=entropy2 >= entropy1, facecolor='#dddddd',
                        interpolate=True)
    ax3[0].set_ylabel('SS entropy local score')
    ax3[0].legend()

    ax3[1].plot(range(N), local_score3, color='black', linewidth=1)
    plt.tight_layout()
    plt.savefig("output/local_score_plot_{}-{}.png".format(pdb_ids[ens1], pdb_ids[ens2]))

def main():
    global M
    pdb_ids, ensembles_features = [], []

    for file in feature_files:
        logging.info("Parsing file: {}".format(file))
        pdb_id = os.path.basename(file)[:-34]
        pdb_ids.append(pdb_id) # get pdb_id from file name
        with open(file) as f:
            data = json.load(f)
            M = len(data.keys())
            logging.info("Computing ensemble features.PDB id: {}".format(pdb_id))
            rg, ss, mrasa, m_rsmd, md, stdev_d = get_ensemble_features(data, pdb_id)
        features = {}
        features['radius_of_giration'] = rg.tolist()
        features['ss_entropy'] = ss.tolist()
        features['median_solvent_accessibility'] = mrasa.tolist()
        features['median_rsmd'] = m_rsmd.tolist()
        features['median_distance'] = md.tolist()
        features['stdev_distance'] = stdev_d.tolist()
        ensembles_features.append(features)
        logging.info("Dumping ensemble features features to file. PDB id: {}".format(pdb_id))
        with open("features/{}_ensemble_features.json".format(pdb_id),'w') as outfile:
            json.dump(features, outfile)

    logging.info("Computing global scores")
    global_score_matrix = get_global_score(pdb_ids, ensembles_features)
    heatmap(pdb_ids, global_score_matrix)

    logging.info("Plotting local score".format(pdb_id))
    for ens1 in range(len(pdb_ids)):
        for ens2 in range(ens1+1,len(pdb_ids)):
            plot_local_score(pdb_ids,ensembles_features,ens1,ens2)
    logging.info("Task 2 program end")


if __name__ == "__main__":
    import argparse, logging, json, os, sys, shutil
    import numpy as np
    from Bio.PDB import Superimposer
    from Bio.PDB.PDBParser import PDBParser
    import matplotlib.pyplot as plt

    # define global variables
    N,M = 0,0

    # configure command line arguments
    parser = argparse.ArgumentParser(description='Code for task 2 Structural Bioinformatics Project')
    parser.add_argument('feature_files', metavar='F', nargs='+',
                        help='Feature files')
    parser.add_argument('--log_stdout', action='store_true', help='Print logging info to standard output')
    parser.add_argument('--reset', action='store_true', help='Empty output, logging and feature folders, except features from task 1')

    # parse arguments
    args = parser.parse_args()
    reset_folders = args.reset
    logging_stdout = args.log_stdout

    # configure logging output
    if logging_stdout:
        logging.basicConfig(stream=sys.stdout, encoding='utf-8', level=logging.INFO,
                            format="%(asctime)s %(levelname)s: %(message)s", filemode='w')
    else:
        if not os.path.exists("logging"):
            os.makedirs("logging")
        logging.basicConfig(filename="logging/task2.log", encoding='utf-8', level=logging.INFO,
                            format="%(asctime)s %(levelname)s: %(message)s", filemode='w')

    # checking that more than one feature file is provided
    if len(args.feature_files) > 1:
        feature_files = args.feature_files

        # checking that required folders exist
        if not os.path.exists("output"): os.makedirs("output")
        if not os.path.exists("features"): os.makedirs("features")

        # cleaning folders if specified
        if reset_folders:
            for folder in ["output", "features"]:
                for filename in os.listdir(folder):
                    file_path = os.path.join(folder, filename)
                    if file_path not in feature_files:
                        try:
                            if os.path.isfile(file_path) or os.path.islink(file_path):
                                os.unlink(file_path)
                            elif os.path.isdir(file_path):
                                shutil.rmtree(file_path)
                        except Exception as e:
                            logging.error('Failed to delete %s. Reason: %s' % (file_path, e))

        logging.info("Program start")
        main()
    else:
        logging.error("Not enough number of feature files (minimum required is 2)")