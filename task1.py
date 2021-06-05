#!/usr/bin/env python

'''0. Utile functions'''
def get_ensemble_dimension(ensemble):
    structure = ensemble[0]
    N,M = 0,len(ensemble)
    for chain in structure:
        for residue in chain:
            if is_aa(residue):  # Filter hetero groups (returns only amino acids)                                      IUPACData.protein_letters_3to1.get(residue.get_resname().capitalize())))
                N += 1
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
            continue
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


def mseDistanceMatrix(matrix1, matrix2):
    array1 = np.asarray(matrix1)
    array2 = np.asarray(matrix2)
    difference_array = np.subtract(array1, array2)
    squared_array = np.square(difference_array)
    mse = squared_array.mean()
    return mse


def meanDifferenceRASA(asa1, asa2):
    array1 = np.asarray(asa1)
    array2 = np.asarray(asa2)
    return np.abs(array1.mean() - array2.mean())


def ssNumberOfDismatch(ss1, ss2):
    str1 = ''.join(ss1)
    str2 = ''.join(ss2)
    val = 1 - SequenceMatcher(None, a=str1, b=str2).ratio()
    return val


def getResidueRGdistanceMatrix():
    dgMat = np.zeros((M, M))
    for i in range(M):
        id = "{}_{}".format(pdb_id, i)
        for j in range(i, M):
            jd = "{}_{}".format(pdb_id, j)
            dgMat[i][j] = np.abs(np.subtract(features[id]['radius_of_giration'], features[jd]['radius_of_giration']))
    u = np.mean(dgMat)
    o = np.std(dgMat)
    return np.divide(dgMat - u, o)


def getResidueASAdistanceMatrix():
    asaMat = np.zeros((M, M))
    for i in range(M):
        id = "{}_{}".format(pdb_id, i)
        for j in range(i, M):
            jd = "{}_{}".format(pdb_id, j)
            asaMat[i][j] = meanDifferenceRASA(features[id]['relative_asas'], features[jd]['relative_asas'])
    u = np.mean(asaMat)
    o = np.std(asaMat)
    return np.divide(asaMat - u, o)


def getResidueSSdistanceMatrix():
    ssMat = np.zeros((M, M))
    for i in range(M):
        id = "{}_{}".format(pdb_id, i)
        for j in range(i, M):
            jd = "{}_{}".format(pdb_id, j)
            ssMat[i][j] = ssNumberOfDismatch(features[id]['secondary_structure'], features[jd]['secondary_structure'])
    u = np.mean(ssMat)
    o = np.std(ssMat)
    return np.divide(ssMat - u, o)


def getResidueDMdistanceMatrix():
    distMat = np.zeros((M, M))
    for i in range(M):
        id = "{}_{}".format(pdb_id, i)
        for j in range(i, M):
            jd = "{}_{}".format(pdb_id, j)
            distMat[i][j] = mseDistanceMatrix(features[id]['distance_matrix'], features[jd]['distance_matrix'])
    u = np.mean(distMat)
    o = np.std(distMat)
    return np.divide(distMat - u, o)

def get_distance_matrix_clustering():
    dg = getResidueRGdistanceMatrix()

    da = getResidueASAdistanceMatrix()

    ds = getResidueSSdistanceMatrix()

    dm = getResidueDMdistanceMatrix()

    return dg, da, ds, dm

def l2_distance(dg, da, ds, dm):
    return np.square(dg ** 2 + da ** 2 + ds ** 2 + dm ** 2)

def l1_distance(dg, da, ds, dm):
    return dg + da + ds + dm

def elbow_clustering(distanceMatrix,elbowErr):
    maxCluster = N-2
    clusters = [[] for x in range(maxCluster)]
    ks = np.zeros(maxCluster)
    err = np.zeros(maxCluster)
    for k in range(2, maxCluster + 2):
        clusters[k - 2], ks[k - 2], err[k - 2] = kmedoids(distance=distanceMatrix, nclusters=int(k), npass=100)
    window_size_smooth=5
    errorSmooth = []
    for i in range(len(ks)):
      frag = ks[np.max([0, i - 1]): np.min([i + window_size_smooth, len(ks)])]
      errorSmooth.append(sum(frag) / len(frag))
    plt.plot(errorSmooth)
    plt.xlabel("Number of clusters")
    plt.ylabel("Smoothed error")
    plt.title("Elbow plot")
    plt.savefig("output/elbow_plot_{}.png".format(pdb_id))
    inertia = []
    inertia_smooth = []
    for i in range(maxCluster - 1):
      inertia.append(np.abs(errorSmooth[i] - errorSmooth[i + 1]))

    window_size_smooth=10
    for i in range(len(inertia)):
      frag = inertia[np.max([0, i - 1]): np.min([i + window_size_smooth, len(inertia)])]
      inertia_smooth.append(sum(frag) / len(frag))

    for i in range(maxCluster - 1):
        if inertia_smooth[i] <= np.mean(inertia_smooth)*elbowErr:
            k = i+2
            break
    cluster=clusters[k-2]
    representatives=[]
    for rep in cluster:
        if rep not in representatives:
            representatives.append(rep)
    return cluster, k, representatives


def graph_heatmap(representatives, heatmap):
    fig, ax = plt.subplots()
    fig.set_figheight(10)
    fig.set_figwidth(10)
    im = ax.imshow(heatmap)

    # We want to show all ticks...
    ax.set_xticks(np.arange(len(representatives)))
    ax.set_yticks(np.arange(len(representatives)))
    # ... and label them with the respective list entries
    ax.set_xticklabels(representatives)
    ax.set_yticklabels(representatives)

    # Rotate the tick labels and set their alignment.
    plt.setp(ax.get_xticklabels(), rotation=45, ha="right",
             rotation_mode="anchor")

    # Loop over data dimensions and create text annotations.
    for i in range(len(representatives)):
        for j in range(len(representatives)):
            text = ax.text(j, i, "%.3f" %(heatmap[i, j]*100),
                           ha="center", va="center", color="w")

    ax.set_title("Graph weight heatmap")
    fig.tight_layout()
    plt.savefig("output/{}_graph_heatmap.png".format(pdb_id))

def graph_printer(nodes,distanceMatrix):
    graph=np.zeros((len(nodes),len(nodes)))
    G = nx.from_numpy_matrix(distanceMatrix)
    labels = nx.get_edge_attributes(G,'weight')
    dge_labels = dict([((u,v,), f"{d['weight']:.2f}") for u,v,d in G.edges(data=True)])
    labels_dict = dict([(u[0], nodes[u[0]]) for u in G.nodes(data=True)])

    plt.figure(figsize=(15,15))
    pos=nx.spring_layout(G)
    nx.draw_networkx_labels(G,pos,labels_dict,font_size=17)
    nx.draw_networkx_nodes(G,pos,nodelist=[j for j in range(0,len(graph))],node_color='lightblue',node_size=1000,ax=None,alpha=0.8)
    nx.draw_networkx_edge_labels(G,pos,edge_labels=dge_labels)
    nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.5)

    nx.write_edgelist(G, "output/text_graph_{}".format(pdb_id))

    plt.savefig('output/{}_graph.png'.format(pdb_id))


def get_variability(list):
    test = np.zeros((len(list), N))
    var = []
    k = 0
    for i in list:
        id = "{}_{}".format(pdb_id, i)
        test[k] = features[id]['relative_asas']
        k = k + 1
    for i in range(N):
        var.append(np.std(test.T[i]))

    var2 = []
    k = 0
    for i in list:
        id = "{}_{}".format(pdb_id, i)
        test[k] = np.sum(features[id]['distance_matrix'], axis=1)
        k = k + 1
    for i in range(N):
        var2.append(np.std(test.T[i]))

    return var


def create_pymol_image(reprentatives):
    out_file = "output/pymol_{}".format(pdb_id)
    cmd.load(pdb_path, pdb_id)
    cmd.split_states(pdb_id, prefix="conf")
    cmd.hide("everything", "all")

    variance = get_variability(reprentatives)
    structures_name = []

    for struct in reprentatives:
        structures_name.append("conf" + ('%04d' % (struct + 1)))
    k = 20
    values = []
    for i in range(0, 100 - k):
        sum = 0
        for i in range(i, i + k):
            sum += variance[i]
        values.append(sum)
    segment_start = np.argmin(values)
    segment_end = segment_start + k

    # string = "show cartoon,"
    struct_string = ""
    for struct in structures_name:
        # string += " " + struct
        #struct_string += "," + struct
    # print(string)
        cmd.show("cartoon", struct)

    strIndx = "& i. {}-{},".format(segment_start, segment_end)
    for struct in structures_name[1:]:
        # str = "align " + struct + " " + strIndx + " " + structures_name[0]
        cmd.align(struct + " " + strIndx, structures_name[0]+ " " + strIndx)
        # print(str)
    # string="zoom "
    cmd.center(structures_name[0]+" "+strIndx)
    cmd.zoom(struct_string)
    # for struct in structures_name:
    # string += " " + struct
    # print("center " + structures_name[0]+ " " +strIndx)
    # print(string)
    # print("png foto" + structures_name[0] +".png")
    features_var = get_variability(range(0, 199))
    norm = colors.Normalize(vmin=np.min(features_var), vmax=np.max(features_var))
    for i, residue in enumerate(Selection.unfold_entities(ensemble[0], "R")):
        rgb = cm.bwr(norm(features_var[i]))
        # print(i, residue.id, structure_rmsd_average[i], rgb)
        cmd.set_color("col_{}".format(i), list(rgb)[:3])
        cmd.color("col_{}".format(i), "resi {}".format(residue.id[1]))
    cmd.png("output/pymol_image_{}.png".format(pdb_id))





def main():
    global pdb_id
    global ensemble
    global N
    global M
    pdb_id = os.path.basename(pdb_path)[:-4]
    ensemble = PDBParser(QUIET=True).get_structure(pdb_id, pdb_path)
    N,M = get_ensemble_dimension(ensemble)
    rg = np.zeros(M)
    global features

    if os.path.isfile("features/{}_single_conformation_features.json".format(pdb_id)):
        logging.info("Loading single conformation features")
        with open("features/{}_single_conformation_features.json".format(pdb_id)) as f:
            features = json.load(f)
    else:
        logging.info("Computing single conformation features")
        for i in range(M):
            id = "{}_{}".format(pdb_id,i)
            features[id] = {}
            logging.info("Structure {} out of {}".format(i,M))
            features[id]['radius_of_giration'] = get_radius_of_gyration(ensemble[i])
            features[id]['relative_asas'] = get_relative_asas(ensemble[i]).tolist()
            features[id]['secondary_structure'] = list(get_secondary_structure(ensemble[i]))
            features[id]['distance_matrix'] =get_distance_matrix(ensemble[i]).tolist()

        logging.info("Dumping single conformation features to file")
        with open("features/{}_single_conformation_features.json".format(pdb_id),'w') as outfile:
            json.dump(features, outfile)

    # logging.info("Showing clustering distance matrix")
    # plt.imshow(distanceMatrix, cmap='hot', interpolation='nearest')
    # plt.show()
    if os.path.isfile("output/clustering_dm_{}.npy".format(pdb_id)):
        logging.info("Loading distance matrix for clustering.")
        distanceMatrix = np.load("output/clustering_dm_{}.npy".format(pdb_id))
    else:
        logging.info("Computing distance matrix for clustering.")
        dg, da, ds, dm = get_distance_matrix_clustering()
        distanceMatrix = l2_distance(dg, da, ds, dm)

        max = np.max(distanceMatrix)
        distanceMatrix = distanceMatrix / max
        for i in range(M):
            logging.info("Structure {} out of {}".format(i, M))
            distanceMatrix[i][i] = 0
            for j in range(i, M):
                distanceMatrix[j][i] = distanceMatrix[i][j]
        #
        if not os.path.exists("output"):
            os.makedirs("output")
        logging.info("Saving distance matrix for clustering.")
        np.save("output/clustering_dm_{}.npy".format(pdb_id), distanceMatrix)

    logging.info("Elbow clustering")
    cluster, k, representatives = elbow_clustering(distanceMatrix, 2)

    repDistanceMatrix = np.zeros((k, k))
    for i in range(k):
        for j in range(k):
            repDistanceMatrix[i][j] = distanceMatrix[representatives[i]][representatives[j]]

    #print(representatives)
    graph_printer(representatives, repDistanceMatrix)
    graph_heatmap(representatives, repDistanceMatrix)

    logging.info("Creating pymol image")



    create_pymol_image(representatives)
    logging.info("Task 1 program end")
    '''
    for i in nodes:
        create_pymol_image(np.where(cluster == i)[0])
    '''
        # cmd.set_color("col_{}".format(i), list(rgb)[:3])
        # cmd.color("col_{}".format(i), "resi {}".format(residue.id[1]))



if __name__ == "__main__":
    import argparse
    import logging
    import json
    import numpy as np
    from Bio.PDB import is_aa, PPBuilder, HSExposureCB
    from Bio.PDB.PDBParser import PDBParser
    from difflib import SequenceMatcher
    from Bio.Cluster import kmedoids
    import math
    import argparse
    import networkx as nx
    from pymol import cmd
    import os
    import sys
    import matplotlib.pyplot as plt
    import warnings
    from Bio.PDB import PDBList, Superimposer, Selection
    from matplotlib import cm, colors
    warnings.filterwarnings("ignore")

    N,M = 0,0
    features = {}
    pdb_id = ""
    ensemble = []


    parser = argparse.ArgumentParser(description='Code for task 1 Structural Bioinformatics Project')
    parser.add_argument('file', metavar='F', nargs=1,
                        help='PDB file path')
    parser.add_argument('--log_stdout', action='store_true', help='Print logging info to standard output')

    args = parser.parse_args()
    pdb_path = args.file[0]
    logging_stdout = args.log_stdout
    if logging_stdout:
        logging.basicConfig(stream=sys.stdout, encoding='utf-8', level=logging.INFO,
                            format="%(asctime)s %(levelname)s: %(message)s", filemode='w')
    else:
        logging.basicConfig(filename="task1.log", encoding='utf-8', level=logging.INFO,
                            format="%(asctime)s %(levelname)s: %(message)s", filemode='w')

    logging.info("Program start")
    main()