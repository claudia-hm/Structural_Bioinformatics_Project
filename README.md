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

### Report 
https://docs.google.com/document/d/1sXajmJBS4c28_oSMq_bLZYVMLse5mOQzaHyH4lh0vFQ/edit?usp=sharing

<!---



!pip install Bio
!apt-get install dssp

### Checklist
This checklist refers to the python files task1.py, task2.p. Only mark as check when definitive working version is available in these python files.

#### Task1
- [x] a. Single conformation features
  - [x] Radius of gyration of the structure
  - [x] Relative accessible surface area (ASA) for each residue.
  - [x] Secondary structure (SS) for each residue
  - [x] Residue distance matrix considering Cα atoms
- [x] b. Graph
- [ ] c. Pymol image

#### Task2
- [x] a. Ensembles features
  - [x] 1. Radius of gyration for each conformation in the ensemble.
  - [x] 2. Secondary structure entropy for each position across ensemble conformations.
  - [x] 3. Median solvent accessibility for each position across ensemble conformations.
  - [x] 4. Median RMSD for each position across ensemble conformations
  - [x] 5. Median distance of each pair of equivalent positions across ensemble conformations.
  - [x] 6. Standard deviation of the distance of each pair of equivalent positions across ensemble
conformations.
- [x] b. Global score heatmap
- [x] c. Plot showing features values along sequence positions

### Commands for running py files

#### Task 1
python task1.py *pdb file*
  
For example:
  
python task1.py data/PED00153e010.pdb --log_stdout

#### Task 2
python task1.py *task1 feature files to compare*
  
For example:
  
python task2.py features/PED00153e007_single_conformation_features.json features/PED00153e008_single_conformation_features.json 
features/PED00153e009_single_conformation_features.json features/PED00153e010_single_conformation_features.json features/PED00153e011_single_conformation_features.json  --log_stdout
-->
