# Structural Bioinformatics Project: Conformational analysis of protein structural ensembles
Claudia Herron Mulet (2029817), Leonardo Amato (2028621) and Matteo Marchesin (1234385)

## Description
This Github repository contains the software required to complete Tasks 1 and 2 of the final project of Structural Bioinformatics.

## Requirements 
The following python modules are required to execute the code:
*argparse, logging, json, os, sys, shutil, warnings, numpy, matplotlib, Bio, math, networkx, difflib, pymol*

## Algorithm description
### Task 1
In this task, we identify the relationships inside a conformational ensemble. First, we calculate the single conformation features as described in the project requirements. Then, to perform clustering and generating a graph, we compute a ```(M,M)``` distance matrix for each of the features, where each cell is the feature distance between structure ```i```and structure ```j```in the ensemble. We combine these 4 matrices into a clustering distance matrix by taking the sum of squares of each feature matrix. Then, we perform clustering using the ```kmedoids```implementation of Biopython and later we identify the optimal number of clusters and hence, the conformational representatives using a technique based on the Elbow plot and inertia. Finally, we compute the pymol image by showing the cluster representatives and aligning them on the residues with less feature variance. Also we color the structures based on the feature variability within the whole ensemble.

### Task 2
In this task 2, we will analyze the relationships between pairs of ensembles. Firs, we calculate the ensemble features as described in the project requirements. Then, we compute the global score applying the formula (4) Lazar et al., 2020 to obtain the *ens_RMS* between pair of ensembles. The output of this section is a global score heatmap between ensembles. Finally, we plot 3 different residue-based local scores: the difference between RSMD, between relative accessible surface area and secondary structure entropy.

More information and discussion available at the [Project Report](https://github.com/claudia-hm/Structural_Bioinformatics_Project/blob/main/SB%20Project%20Report.pdf).

## Warning 
These programs output the feature files inside the ```features``` folder and other required and complementary results in the ```output```. The code is not very fast. All features and matrices have been precomputed and they are automatically loaded if not existent in the ```features``` and  ```output``` folders. To visualize in real time the computations that are being performed is recommended the use of the flag ```--log_stdout``` (later explained). It is possible to force the recomputing of the files using the ```--reset``` flag.

## Usage 
1. Download this repository 
2. ```cd  Structural_Bioinformatics_Project```
3. Load ensemble pdb files into ```data``` directory. The files must conserve the original name from [PED](https://proteinensemble.org), e.g., ```PED00153e007.pdb```
4. Run ```task1.py``` providing as input the path to the PDB file: ```python task1.py pdb_file```, e.g,  ```python task1.py data/PED00153e007.pdb```
5. Run ```task2.py``` providing as input the path to the feature files generated in Task 1: ```python task2.py feature_file_1 feature_file_2 ...```, e.g., ```python task2.py features/PED00153e007_single_conformation_features.json features/PED00153e008_single_conformation_features.json```
6. Both python files output their results to the folders ```output```and ```features```folders. In addition, they can be called using the following flags:
* ```--log_stdout```: redirects the logging information to the standard output. When not present, the logging information is printed to the files ```task1.log```and ```task2.log```.
* ```--reset```: when present, the program cleans the content of the ```features``` and ```output``` folders, except those needed files for the execution.
* ```-h```: help


