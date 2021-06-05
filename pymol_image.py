import numpy as np
from Bio.PDB import PDBList, Superimposer, Selection
from Bio.PDB.PDBParser import PDBParser
import pymol
from pymol import cmd
from matplotlib import cm, colors
import os



pdb_id = "PED00153e007"
pdb_file = "data/PED00153e007.pdb"
ensemble = PDBParser(QUIET=True).get_structure(pdb_id, pdb_file)


out_file = "output/pymol_{}".format(pdb_id)

cmd.load(pdb_file, pdb_id)

#
cmd.split_states(pdb_id, prefix="conf")
#cmd.delete(pdb_id)
#
cmd.hide("everything","all")
cmd.show("cartoon", "conf0014")
cmd.show("cartoon", "conf0072")
cmd.show("cartoon", "conf0150")
#cmd.hide(pdb_id)  # Hide lines
cmd.center("conf0014")
cmd.zoom("conf0014")

features_var = [0.2235496810038598, 0.21423024627626014, 0.2581930113583977, 0.24264139511880942, 0.200133399455409, 0.22433974815622743, 0.21810326172969113, 0.1783636551783323, 0.2113018718292807, 0.22670356288721857, 0.22669942616517647, 0.24442025488267605, 0.22161663097867193, 0.20037311208054298, 0.20355431180277345, 0.21932138496012332, 0.20429591596944602, 0.2053923025290539, 0.21936380155094298, 0.2290554018522377, 0.22070719646343404, 0.2103757677444239, 0.22007199369091793, 0.26027961238509106, 0.26468787907766106, 0.24227121066766014, 0.25567827446291647, 0.24757045693067103, 0.2678658275010321, 0.25824813416785536, 0.24216752348590187, 0.2721757333607947, 0.23051065833970877, 0.25170478270042623, 0.2408347665782758, 0.19439831424223555, 0.20526215489648247, 0.20701017538452338, 0.16375526612461055, 0.19099841608736012, 0.17713628662953287, 0.20331905052743574, 0.1971413152512803, 0.1971769560869296, 0.20886334410385443, 0.23570745620112005, 0.23148967568485912, 0.23155989671260505, 0.24246047783275354, 0.2688215527488618, 0.259073908419263, 0.2268229048443392, 0.2830922793042832, 0.23751948508406678, 0.2822854204543102, 0.2552737035167884, 0.25064584843077414, 0.25996579553306454, 0.23584158186043772, 0.25525646230087856, 0.24716008479276966, 0.24419044860383587, 0.23196153431609234, 0.23267924841749937, 0.18912028241936132, 0.24031905462232553, 0.20145095772786295, 0.20759493741292664, 0.23592334092418854, 0.2334701028354617, 0.2714229056638644, 0.22578503719481952, 0.248419981910497, 0.2498178348601916, 0.23498385180053077, 0.2646478283151562, 0.25247466985995237, 0.24634473576910426, 0.2628718672824459, 0.24271608693837549, 0.26061101868073816, 0.2086088193186364, 0.24411251728781796, 0.2290760758386054, 0.2545907084572669, 0.2759324962916536, 0.258307262893001, 0.2578367017010391, 0.25292120144679003, 0.23525033378358207, 0.27236462156511254, 0.22368279205132335, 0.2175750961533168, 0.22262839212927368, 0.2401190377089321, 0.2103311802417775, 0.1975805159981368, 0.19358278420632355, 0.1997621311567337, 0.16885900734998077]
argmin_var = np.argmin(features_var)
cmd.align("conf0072 & i. 38", "conf0014 & i. 38")
cmd.align("conf0150 & i. 38", "conf0014 & i. 38")

norm = colors.Normalize(vmin=min(features_var), vmax=max(features_var))
for i, residue in enumerate(Selection.unfold_entities(ensemble[14], "R")):
    rgb = cm.bwr(norm(features_var[i]))
    # print(i, residue.id, structure_rmsd_average[i], rgb)
    cmd.set_color("col_{}".format(i), list(rgb)[:3])
    cmd.color("col_{}".format(i), "resi {}".format(residue.id[1]))

cmd.png(out_file, width=2000, height=2000, ray=1)






# pdb_id = '2k0e'
# pdb_file = "data/pdb{}.ent".format(pdb_id)
# out_file = "data/pymol_{}".format(pdb_id)
# window_size = 9
#
#
# # Fetch a PDB file to the current dir
# pdbl = PDBList()
# pdbl.retrieve_pdb_file(pdb_id, pdir='data/', file_format='pdb') # Will save to pdbXXXX.ent
#
# # Load the structure
# structure = PDBParser(QUIET=True).get_structure(pdb_id, pdb_file)
#
# # Superimpose all models to the first model, fragment-by-fragment (sliding window)
# super_imposer = Superimposer()
# structure_rmsd_fragments = []  # RMSD, no_models X no_fragments X fragment_size
# ref_model = [atom for atom in structure[0].get_atoms() if atom.get_name() == "CA"]  # CA of the first model
#
# # Iterate all models
# for i, model in enumerate(structure):
#     if i > 0:
#         model_rmsd = []  # RMSD, no_fragment X fragment_size
#         alt_model = [atom for atom in model.get_atoms() if atom.get_name() == "CA"]  # coords of the model
#
#         # Iterate fragments
#         for start in range(len(ref_model) - window_size):
#             end = start + window_size
#             ref_fragment = ref_model[start:end]
#             alt_fragment = alt_model[start:end]
#
#             # Calculate rotation/translation matrices
#             super_imposer.set_atoms(ref_fragment, alt_fragment)
#             # print(super_imposer.rms, super_imposer.rotran)
#
#             # Rotate-translate coordinates
#             alt_fragment_coord = np.array([atom.get_coord() for atom in alt_fragment])
#             alt_fragment_coord = np.dot(super_imposer.rotran[0].T, alt_fragment_coord.T).T
#             alt_fragment_coord = alt_fragment_coord + super_imposer.rotran[1]
#
#             # Calculate RMSD
#             # https://en.wikipedia.org/wiki/Root-mean-square_deviation_of_atomic_positions
#             ref_fragment_coord = np.array([atom.get_coord() for atom in ref_fragment])
#             dist = ref_fragment_coord - alt_fragment_coord
#             # rmsd_fragment = np.sqrt(np.sum(dist * dist) / window_size)  # Total RMSD of the fragment. Identical to super_imposer.rms
#             rmsd_res = np.sqrt(np.sum(dist * dist, axis=1))  # RMSD for each residue of the fragment
#
#             model_rmsd.append(rmsd_res)
#
#         structure_rmsd_fragments.append(model_rmsd)
#
#
# # Calculate average RMSD per position
# structure_rmsd_fragments = np.array(structure_rmsd_fragments)  # no_models X no_fragments X fragment_size
# # Calculate the RMSD average for each fragments along all models
# structure_rmsd_fragments = np.average(structure_rmsd_fragments, axis=0)  # no_fragments X fragment_size
# # Pad with right zeros to reach the sequence length (no_fragments + fragment_size)
# structure_rmsd_fragments = np.pad(structure_rmsd_fragments, ((0, 0), (0, structure_rmsd_fragments.shape[0])))
# print(structure_rmsd_fragments.shape, len(ref_model))
#
# # Roll the fragments one by one (add heading zeros)
# for i, row in enumerate(structure_rmsd_fragments):
#     structure_rmsd_fragments[i] = np.roll(row, i)
#
#
# # Calculate average along columns of overlapping fragments (average RMSD per residue)
# structure_rmsd_average = np.average(structure_rmsd_fragments, axis=0)
# print(structure_rmsd_average.shape)
#
#
# #################### PYMOL #####################
#
# #pymol.finish_launching()  # Open Pymol
#
# cmd.load(pdb_file, pdb_id)  # Download the PDB
# # cmd.load("data/pdb{}.ent".format(pdb_id), pdb_id)  # Load from file
#
# cmd.remove("resn hoh")  # Remove water molecules
# cmd.hide("lines", "all")  # Hide lines
# cmd.show("cartoon", pdb_id)  # Show cartoon
#
# norm = colors.Normalize(vmin=min(structure_rmsd_average), vmax=max(structure_rmsd_average))
# for i, residue in enumerate(Selection.unfold_entities(structure[0], "R")):
#     rgb = cm.bwr(norm(structure_rmsd_average[i]))
#     # print(i, residue.id, structure_rmsd_average[i], rgb)
#     cmd.set_color("col_{}".format(i), list(rgb)[:3])
#     cmd.color("col_{}".format(i), "resi {}".format(residue.id[1]))
#
# cmd.png(out_file, width=2000, height=2000, ray=1)
