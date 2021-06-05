#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun  3 14:24:04 2021

@author: matteo
"""









import numpy as np
#from Bio.PDB import PDBList, Superimposer, is_aa
#from Bio.PDB.PDBParser import PDBParser
import requests
import math
#import matplotlib.pyplot as plt
import pymol
from pymol import cgo, cmd, util
import os

os.getcwd()

pdb_id = ''

#cmd.fetch(pdb_id, pdb_id, path="data/")  # Download the PDB
# cmd.load("data/pdb{}.ent".format(pdb_id), pdb_id)  # Load from file
ped_id = "PED00153"
url = "https://proteinensemble.org/api/" + ped_id
resp_json = requests.get(url).json()
print(resp_json["title"])
ensembles_ids = []
for curr_ensemble in resp_json["ensembles"]:
    ensembles_ids.append(curr_ensemble["ensemble_id"])
    
    
    
pdb_dir ="data/"
pdb_file ="{}.pdb".format(ensembles_ids[0])
pdb_id = ensembles_ids[0]

pymol.finish_launching()  # Open Pymol

cmd.load(pdb_dir+pdb_file, pdb_id) 
#cmd.load(pdb_dir+"57_"+pdb_file, pdb_id+"_57") 

cmd.remove("resn hoh")  # Remove water molecules
cmd.hide("lines", "all")  # Hide lines
cmd.show("cartoon", pdb_id)  # Show cartoon
cmd.hide("sticks", "hetatm")  # Show hetero atoms as sticks
cmd.hide()
cmd.split_states(pdb_id)



cmd.show('cartoon',pdb_id+'_0172')
cmd.show('cartoon',pdb_id+'_0002')
cmd.set_color(pdb_id+'_0002', 'red')


atom1= 'A/58/CA'
#prendo le coordinate dell'atomo 1 definito sopra per il model 172
atm1_172=cmd.get_coords(atom1,172)

atm1_2=cmd.get_coords(atom1,2)
cmd.align(pdb_id+'_0172',pdb_id+'_0002')


