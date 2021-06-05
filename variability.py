import json
import numpy as np

pdb_id = "PED00153e007"
file = "features/PED00153e007_single_conformation_features.json"
with open(file) as f:
    features = json.load(f)
    M = len(features.keys())
N = 100

### pymol part ###
rasas = []
for i in range(M):
    id = "{}_{}".format(pdb_id, i)
    rasas.append(np.array(features[id]['relative_asas']))
rasas = np.array(rasas)
# print(rasas.shape)
# print(len(rasas[:,0]))
rasas_std = np.zeros(N)
for i in range(N):
     rasas_std[i] = np.std(rasas[:,i])


