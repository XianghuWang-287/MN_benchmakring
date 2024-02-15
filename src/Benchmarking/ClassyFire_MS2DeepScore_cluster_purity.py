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
import random

# classifier_df = pd.read_csv("./save_rosetta_ClassyfireOverhaul.csv",sep=',')
# library = "GNPS-SELLECKCHEM-FDA-PART2"
# summary_file_path = "./data/summary/"+library+"_summary.tsv"
# merged_pairs_file_path = "./data/merged_paris/"+library+"_merged_pairs.tsv"
# cluster_summary_df = pd.read_csv(summary_file_path)
# all_pairs_df = pd.read_csv(merged_pairs_file_path, sep='\t')

# for index, row in cluster_summary_df.iterrows():
#     print("Spectrum ID:", row["spectrum_id"])
#     if(classifier_df["Lib_ids"].isin([row["spectrum_id"]]).any()):
#         print("library smiles:",row['Smiles'],"classifier df smiles:",classifier_df.loc[classifier_df["Lib_ids"]==row["spectrum_id"]]["SMILES"].iloc[0])
#     else:
#         print("missing data")

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
def calculate_correct_classification_percentage(components, node_types_df, threshold=0.7):
    correctly_classified_clusters = 0
    total_clusters = len(components)

    for component in components:
        cluster_nodes = list(component.nodes())
        node_index = [node - 1 for node in cluster_nodes]
        try:
            cluster_types = node_types_df.loc[node_index, 'Superclass']
            cluster_types = cluster_types.replace({np.nan: "Unknown"})
            dominant_type = cluster_types.value_counts().idxmax()
            if dominant_type == "Unknown":
                continue
            percentage = cluster_types.value_counts()[dominant_type] / len(cluster_nodes)

            if percentage >= threshold:
                correctly_classified_clusters += 1
        except Exception as e:
            continue
            #print(e)
            #print(node_index)

    return correctly_classified_clusters / total_clusters

def edge_correct_classification_percentage(components,node_types_df):
    correct_num = 0
    total_number = 0
    for component in components:
        total_number = total_number + component.number_of_edges()
        for edge in component.edges():
            node1 = edge[0]
            node2 = edge[1]
            node1_type = node_types_df.iloc[node1-1]['Kingdom']
            node2_type = node_types_df.iloc[node2-1]['Kingdom']
            if pd.isna(node1_type) or pd.isna(node2_type):
                continue
            if node1_type == node2_type:
                correct_num = correct_num+1
    return correct_num/total_number

def CalN50(lst,toatl_num, percentage):
    # Sort the list in descending order
    sorted_list = sorted(lst, reverse=True)



    # Calculate the target sum (20% of the total sum)
    target_sum = percentage * toatl_num

    # Initialize variables
    current_sum = 0
    min_elements = sorted_list[0]

    # Iterate through the sorted list and add elements until the target sum is reached
    for num in sorted_list:
        current_sum += num
        min_elements = num
        if current_sum >= target_sum:
            break

    return min_elements


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Using realignment method to reconstruct the network')
    parser.add_argument('-m', type=str, required=True, default="spec.mgf", help='input mgf filename')
    parser.add_argument('--input', type=str,required=True,default="input_library", help='input libraries')
    args = parser.parse_args()
    input_lib_file = args.input
    merged_file = args.m
    summary_file_path = "./data/summary/"+input_lib_file.replace(".mgf", "") + "_summary.tsv"
    classifer_file_path = "./" + input_lib_file.replace(".mgf", "") + "_classifier_results.csv"
    cluster_summary_df = pd.read_csv(summary_file_path)
    classifer_df = pd.read_csv(classifer_file_path)
    merged_file = merged_file+input_lib_file
    threshold_list = [0.9,0.91,0.92,0.93,0.94,0.95,0.96,0.97,0.98,0.99,0.991,0.992,0.993,0.994,0.995,0.996,0.997,0.998,0.999,0.9995,0.9996,0.9997,0.9998,0.9999]
    N20_list=[]
    score_list=[]
    merged_file_original = "./data/merged_paris/"+input_lib_file.replace(".mgf", "") + "_merged_pairs.tsv"
    original_all_pairs_df = pd.read_csv(merged_file_original, sep='\t')
    G_original = nx.from_pandas_edgelist(original_all_pairs_df, "CLUSTERID1", "CLUSTERID2", "Cosine")
    num_of_nodes = G_original.number_of_nodes()
    print(num_of_nodes)
    purity_score_list = []
    N50_list = []
    for threshold in threshold_list:
        merged_pairs_file_path = convert_extension(merged_file,threshold,"tsv")
        print(merged_pairs_file_path)
        all_pairs_df = pd.read_csv(merged_pairs_file_path, sep='\t')
        G_all_pairs = nx.from_pandas_edgelist(all_pairs_df, "CLUSTERID1", "CLUSTERID2", "Cosine")
        print('graph with {} nodes and {} edges'.format(G_all_pairs.number_of_nodes(), G_all_pairs.number_of_edges()))
        print("constructing dic for finger print")
        components_all_list = []
        score_all_pairs_filter_list=[]
        components = [G_all_pairs.subgraph(c).copy() for c in nx.connected_components(G_all_pairs)]
        correct_classification_percentage = calculate_correct_classification_percentage(components, classifer_df)
        purity_score_list.append(correct_classification_percentage)
        components_number = [x.number_of_nodes() for x in components]
        N50_list.append(CalN50(components_number,G_all_pairs.number_of_nodes(), 0.2))
    print(purity_score_list)
    print(N50_list)