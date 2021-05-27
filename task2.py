#!/usr/bin/env python

def main():

    ## not working
    '''
    for file in feature_files:
        with open(file) as f:
            data = json.load(f)
            print(data.keys())
    '''





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
