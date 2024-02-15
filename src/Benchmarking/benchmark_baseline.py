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
    number=dfnew.values[0][1]
    row_old = dfnew.values[0][1]
    if len(dfnew.values) ==1:
        return row_old
    for row in dfnew.values[1:]:
        if (number >= node_numbers*N_ratio):
            return row_old
        else:
            number=number+ row[1]
            row_old = row[1]
    return 1

def convert_extension(input_string,threshold, new_extension):
    # Split the input string into name and extension
    name, old_extension = input_string.rsplit('.', 1)
    name = name + "_" + str(threshold)

    # Check if the old extension is "mgf"
    if old_extension.lower() == 'mgf':
        # Concatenate the name with the new extension
        output_string = f"{name}.{new_extension}"
        return output_string
    else:
        # If the input string doesn't end with ".mgf", return an error message or handle it as needed
        return "Invalid input file format"
# def cal_N50(df, node_numbers,N_ratio):
#     dfnew=df.sort_values('number',ascending=False)
#     number=0
#     row_old = dfnew.values[0][1]
#     for row in dfnew.values:
#         if (number >= node_numbers*N_ratio):
#             return row_old
#         else:
#             number=number+ row[1]
#             row_old = row[1]
#     return 1
if __name__ == '__main__':
    #pass arguments
    parser = argparse.ArgumentParser(description='Using realignment method to reconstruct the network')
    parser.add_argument('-m', type=str, required=True, default="spec.mgf", help='input mgf filename')
    parser.add_argument('--input', type=str,required=True,default="input_library", help='input libraries')
    args = parser.parse_args()
    input_lib_file = args.input
    merged_file = args.m

    #read libraries from input file
    # with open(input_lib_file,'r') as f:
    #     libraries = f.readlines()

    summary_file_path = "./data/summary/"+input_lib_file.replace(".mgf", "") + "_summary.tsv"
    # merged_pairs_file_path = "./data/merged_paris/"+library+"_merged_pairs.tsv"
    cluster_summary_df = pd.read_csv(summary_file_path)
    merged_file = merged_file+input_lib_file
    print(merged_file)
    threshold_list = [0.5,0.6,0.7,0.8,0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.991,0.992,0.993,0.994,0.995,0.996,0.997]
    N20_list=[]
    score_list=[]
    merged_file_original = "./data/merged_paris/"+input_lib_file.replace(".mgf", "") + "_merged_pairs.tsv"
    original_all_pairs_df = pd.read_csv(merged_file_original, sep='\t')
    G_original = nx.from_pandas_edgelist(original_all_pairs_df, "CLUSTERID1", "CLUSTERID2", "Cosine")
    num_of_nodes = G_original.number_of_nodes()
    print(num_of_nodes)
    for threshold in threshold_list:
        merged_pairs_file_path = convert_extension(merged_file,threshold,"tsv")
        print(merged_pairs_file_path)
        all_pairs_df = pd.read_csv(merged_pairs_file_path, sep='\t')
        G_all_pairs = nx.from_pandas_edgelist(all_pairs_df, "CLUSTERID1", "CLUSTERID2", "Cosine")
        print('graph with {} nodes and {} edges'.format(G_all_pairs.number_of_nodes(), G_all_pairs.number_of_edges()))
        print("constructing dic for finger print")
        dic_fp = fingerprint_dic_construct(cluster_summary_df)
        results_df_list=[]
        score_all_pairs_filter_list=[]
        components = [G_all_pairs.subgraph(c).copy() for c in nx.connected_components(G_all_pairs)]
        for component in tqdm(components):
            score_all_pairs_filter_list.append(subgraph_score_dic(component,cluster_summary_df,dic_fp))
        all_pairs_filter_number = [len(x) for x in components]
        df_all_pairs_filter = pd.DataFrame(list(zip(score_all_pairs_filter_list, all_pairs_filter_number)),columns=['score', 'number'])
        results_df_list.append(df_all_pairs_filter)
        # result_file_path = "./results-base/"+library+"_baseline_benchmark.pkl"
        # with open(result_file_path, 'wb') as file:
        #     pickle.dump(results_df_list, file)

        print(np.array([cal_N50(x, num_of_nodes, 0.2) for x in results_df_list]))
        print(np.array([weighted_average(x, 'score', 'number') for x in results_df_list]))
        N20_list.append(np.array([cal_N50(x, num_of_nodes, 0.2) for x in results_df_list])[0])
        score_list.append(np.array([weighted_average(x, 'score', 'number') for x in results_df_list])[0])
    print(N20_list)
    print(score_list)




