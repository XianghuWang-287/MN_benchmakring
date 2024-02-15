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
def cal_N50(df, node_numbers,N_ratio):
    dfnew=df.sort_values('number',ascending=False)
    number=0
    for row in dfnew.values:
        if (number >= node_numbers*N_ratio):
            return row_old
        else:
            number=number+ row[1]
            row_old = row[1]


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
        merged_pairs_file_path = "./data/merged_paris/"+library+"_merged_pairs.tsv"
        cluster_summary_df = pd.read_csv(summary_file_path)
        all_pairs_df = pd.read_csv(merged_pairs_file_path, sep='\t')
        G_all_pairs = nx.from_pandas_edgelist(all_pairs_df, "CLUSTERID1", "CLUSTERID2", "Cosine")
        print('graph with {} nodes and {} edges'.format(G_all_pairs.number_of_nodes(), G_all_pairs.number_of_edges()))
        print("constructing dic for finger print")
        dic_fp = fingerprint_dic_construct(cluster_summary_df)
        x_max_number = [x for x in range(1, 40, 2)]
        y_max_number = [y for y in range(2, 402, 20)]
        y_weight_avg = []
        components_all_list = []
        results_df_list = []
        G_all_pairs_filter = G_all_pairs.copy()
        filter_top_k(G_all_pairs_filter, 10)
        filter_component(G_all_pairs_filter, 100)
        G_all_pairs_filter_copy = G_all_pairs_filter.copy()
        for node in G_all_pairs_filter_copy.nodes():
            if G_all_pairs_filter.degree[node] == 0:
                G_all_pairs_filter.remove_node(node)
        score_all_pairs_filter_list = []
        components = [G_all_pairs_filter.subgraph(c).copy() for c in nx.connected_components(G_all_pairs_filter)]
        components_all_list.append(components)
        for component in tqdm(components):
            score_all_pairs_filter_list.append(subgraph_score_dic(component,cluster_summary_df,dic_fp))
        all_pairs_filter_number = [len(x) for x in components]
        df_all_pairs_filter = pd.DataFrame(list(zip(score_all_pairs_filter_list, all_pairs_filter_number)),columns=['score', 'number'])
        results_df_list.append(df_all_pairs_filter)
        print(np.array([cal_N50(x, 1093, 0.2) for x in results_df_list]))
        print(np.array([weighted_average(x, 'score', 'number') for x in results_df_list]))





