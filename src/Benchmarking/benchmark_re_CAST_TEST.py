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
def all_elements_duplicated(list_of_lists):
    flattened_list = [item for sublist in list_of_lists for item in sublist]

    unique_set = set(flattened_list)

    if len(flattened_list) == len(unique_set):
        return False  # Not all elements are duplicates
    else:
        return True  # All elements are duplicates

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
        re_align_edge_file_path = "./alignment_results/"+library+"_realignment.pkl"
        cluster_summary_df = pd.read_csv(summary_file_path)
        all_pairs_df = pd.read_csv(merged_pairs_file_path, sep='\t')
        G_all_pairs = nx.from_pandas_edgelist(all_pairs_df, "CLUSTERID1", "CLUSTERID2", "Cosine")
        print('graph with {} nodes and {} edges'.format(G_all_pairs.number_of_nodes(), G_all_pairs.number_of_edges()))
        print("constructing dic for finger print")
        dic_fp = fingerprint_dic_construct(cluster_summary_df)
        print("constructing the re-alignment graph")
        G_all_pairs_realignment = G_all_pairs.copy()
        with open(re_align_edge_file_path, 'rb') as f:
            realignment_edgelist = pickle.load(f)
        new_list=[]
        for item in realignment_edgelist:
            if(item != None):
                if item[2]>=0.8:
                    new_list.append(item)
        for item in new_list:
            if (item != None):
                G_all_pairs_realignment.add_edge(item[0], item[1], Cosine=item[2])
        results_df_list = []
        thresholds = [x / 100 for x in range(75, 95)]
        for threshold in thresholds:
            cast_cluster = CAST_cluster(G_all_pairs_realignment, threshold)
            cast_score_list = []
            cast_components = [G_all_pairs_realignment.subgraph(c).copy() for c in cast_cluster]
            print(all_elements_duplicated(cast_cluster))
            re_align_set=set(cast_components[0].edges())
            for component in cast_components:
                re_align_set = re_align_set.union(set(component.edges()))
            original_set = set(G_all_pairs.edges())
            common_edges =  re_align_set.intersection(original_set)
            new_add = re_align_set - common_edges
            print("re_align_set:{}".format(len(re_align_set)))
            print("originial_set:{}".format(len(original_set)))
            print("common_edges:{}".format(len(common_edges)))
            G_intersection = nx.Graph()
            G_intersection.add_edges_from(common_edges)
            # components_num = nx.number_connected_components(G_intersection)
            # components = [G_intersection.subgraph(c).copy() for c in nx.connected_components(G_intersection)]
            # plt.figure(figsize=(18, 100))
            # for component in components:
            #     plt.subplot(components_num // 8 + 1, 8, components.index(component) + 1)
            #     nx.draw(component, node_size=2)
            # plt.show()
            score_intersection_list=[]
            components = [G_intersection.subgraph(c).copy() for c in nx.connected_components(G_intersection)]
            for component in tqdm(components):
                score_intersection_list.append(subgraph_score_dic(component, cluster_summary_df, dic_fp))
            intersection_number = [len(x) for x in components]
            df_intersection = pd.DataFrame(list(zip(score_intersection_list, intersection_number)),
                                               columns=['score', 'number'])
            results_df_list.append(df_intersection)
        result_file_path = "./results-intersection/" + library + "_intersection.pkl"
        with open(result_file_path, 'wb') as file:
            pickle.dump(results_df_list, file)












