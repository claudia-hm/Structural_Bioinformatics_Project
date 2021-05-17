# Project understanding
A markdown file where we can explain related concepts to the project and hints.

## Glossary
* **Intrinsically disordered protein (IDP)**: a protein that lacks a fixed or ordered three-dimensional structure, typically in the absence of its macromolecular interaction partners, such as other proteins or RNA.
* **Conformational ensembles (aka Structural ensembles)**: experimentally constrained computational models describing the structure of intrinsically unstructured proteins.
* **Conformations**: the different 3-dimensional arrangements that the molecule can acquire by freely rotating around σ-bonds.
* **Sigma bonds (σ bonds)**: the strongest type of covalent chemical bond.
* **Solving IDP/IDR ensembles**: finding a set of IDP conformations that is a is a faithful representation of the real physical state of the IDP/IDR.
* **Radius of gyration**: the root-mean-square average of the distance of all scattering elements from the center of mass of the molecule. It indicates the spread of the coordinates of the protein.
* **Solvent exposure**: of an amino acid in a protein measures to what extent the amino acid is accessible to the solvent (usually water) surrounding the protein
* **Relative accessible surface area or relative solvent accessibility**: a measure of residue solvent exposure.

## Directly related lectures
* **Lecture 15**: 
  1. Explanation of project requirements. Read articles. Idea of the project: pairwise comparison between the different ensebles in PED. We should provide biological considerations about the ensemble in the report (second article!). Visualize one ensemble in pymol and let the user understand what is happening.
  2. Code:
      *  secondary_structure.py: explanation of dssp. Dssp outputs ***relative ASA*** (if it is higher than a threshold then it is exposed). There is the code for calculating ***phi-psi angles***.
      * radius_gyration.py: TO-DO
      * superposition_fragments.py: TO-DO
  3. Introduction to IDPs and NMR.
* **Lecture 16**: 
  1. TM score and TM align: not needed for the project.
  2. Code:
      * pymol_states.py: TO-DO
      * pymol_lines.py: TO-DO
* **Lecture 17**: theory related to IDP.
* **Lecture 18**: theory related to IDP.
