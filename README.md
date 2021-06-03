# Structural_Bioinformatics_Project

### Report 
https://docs.google.com/document/d/1sXajmJBS4c28_oSMq_bLZYVMLse5mOQzaHyH4lh0vFQ/edit?usp=sharing


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
- [ ] b. Graph
- [ ] c. Pymol image

#### Task2
- [x] a. Ensembles features
  - [x] 1. Radius of gyration for each conformation in the ensemble.
  - [ ] 2. Secondary structure entropy for each position across ensemble conformations.
  - [x] 3. Median solvent accessibility for each position across ensemble conformations.
  - [ ] 4. Median RMSD for each position across ensemble conformations
  - [x] 5. Median distance of each pair of equivalent positions across ensemble conformations.
  - [x] 6. Standard deviation of the distance of each pair of equivalent positions across ensemble
conformations.
- [x] b. Global score heatmap
- [ ] c. Plot showing features values along sequence positions

### Commands for running py files

#### Task 1
python task1.py *pdb file*
  
For example:
  
python task1.py data/PED00153e010.pdb

#### Task 2
python task1.py *task1 feature files to compare*
  
For example:
  
python task2.py PED00153e007_single_conformation_features.json PED00153e008_single_conformation_features.json PED00153e009_single_conformation_features.json PED00153e010_single_conformation_features.json PED00153e011_single_conformation_features.json
