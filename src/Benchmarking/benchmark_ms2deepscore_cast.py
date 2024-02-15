import pandas as pd
import networkx as nx
import matplotlib.pyplot as plt
import time
import urllib
from tqdm import tqdm
import matplotlib.colors as colors
import matplotlib.cbook as cbook
from matplotlib import cm
import matplotlib as mpl
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Descriptors import ExactMolWt
from rdkit.Chem.Draw import MolToFile
from rdkit.DataStructs import FingerprintSimilarity
from rdkit.Chem.Fingerprints.FingerprintMols import FingerprintMol
from rdkit.Chem.rdMolDescriptors import CalcMolFormula
import requests
import plotly.express as px
from IPython.display import HTML
from utility import *
from pyteomics import mgf
from Classic import *
from multiprocessing import Pool
import collections
from typing import List, Tuple
import pickle
import os
import argparse



if __name__ == '__main__':
    #pass arguments
    parser = argparse.ArgumentParser(description='Using realignment method to reconstruct the network')
    parser.add_argument('--input', type=str,required=True,default="input_library.txt", help='input libraries')
    args = parser.parse_args()
    input_lib_file = args.input

    #read libraries from input file
    with open(input_lib_file,'r') as f:
        libraries = f.readlines()

    for library in libraries:
        library = library.strip('\n')
        print("starting benchmarking library:"+library)
        summary_file_path = "./data/summary/"+library+"_summary.tsv"
        merged_pairs_file_path = "./data/raw_ms2deepscore/"+library+"_0.9.tsv"
        cluster_summary_df = pd.read_csv(summary_file_path)
        all_pairs_df = pd.read_csv(merged_pairs_file_path, sep='\t')
        G_all_pairs = nx.from_pandas_edgelist(all_pairs_df, "CLUSTERID1", "CLUSTERID2", "Cosine")
        print('graph with {} nodes and {} edges'.format(G_all_pairs.number_of_nodes(), G_all_pairs.number_of_edges()))
        print("constructing dic for finger print")
        dic_fp = fingerprint_dic_construct_InCHI(cluster_summary_df)
        y_weight_avg = []
        components_all_list = []
        results_df_list = []
        thresholds = [x / 1000 for x in range(900, 999,5)]
        for threshold in tqdm(thresholds):
            cast_cluster = CAST_cluster(G_all_pairs, threshold)
            cast_score_list = []
            cast_components = [G_all_pairs.subgraph(c).copy() for c in cast_cluster]
            for component in cast_components:
                cast_score_list.append(subgraph_score_dic(component, cluster_summary_df, dic_fp))
            cast_number = [len(x) for x in cast_cluster]
            df_cast = pd.DataFrame(list(zip(cast_score_list, cast_number)), columns=['score', 'number'])
            results_df_list.append(df_cast)
        result_file_path = "./results_cast_ms2deepscore/" + library + "_cast_benchmark.pkl"
        with open(result_file_path, 'wb') as file:
            pickle.dump(results_df_list, file)




