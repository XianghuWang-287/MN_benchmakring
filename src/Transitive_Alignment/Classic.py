import pandas as pd
# from gnpsdata import taskresult
# from gnpsdata import workflow_classicnetworking
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


def prune_component(G, component, cosine_delta=0.02):
    component_edges = get_edges_of_component(G, component)

    min_score = 1000
    for edge in component_edges:
        min_score = min(min_score, edge[2]["Cosine"])

    cosine_threshold = cosine_delta + min_score
    for edge in component_edges:
        if edge[2]["Cosine"] < cosine_threshold:
            #print(edge)
            G.remove_edge(edge[0], edge[1])
def get_edges_of_component(G, component):
    component_edges = {}
    for node in component:
        node_edges = G.edges((node), data=True)
        for edge in node_edges:
            if edge[0] < edge[1]:
                key = str(edge[0]) + "-" + str(edge[1])
            else:
                key = str(edge[1]) + "-" + str(edge[0])
            component_edges[key] = edge

    component_edges = component_edges.values()
    return component_edges
def filter_component(G, max_component_size):
    if max_component_size == 0:
        return

    big_components_present = True

    while big_components_present == True:
        big_components_present = False
        components = nx.connected_components(G)
        for component in components:
            if len(component) > max_component_size:
                prune_component(G, component)
                big_components_present = True
        #print("After Round of Component Pruning", len(G.edges()))
def filter_top_k(G, top_k):
    #print("Filter Top_K", top_k)
    #Keeping only the top K scoring edges per node
    #print("Starting Numer of Edges", len(G.edges()))

    node_cutoff_score = {}
    for node in G.nodes():
        node_edges = G.edges((node), data=True)
        node_edges = sorted(node_edges, key=lambda edge: edge[2]["Cosine"], reverse=True)

        edges_to_delete = node_edges[top_k:]
        edges_to_keep = node_edges[:top_k]

        if len(edges_to_keep) == 0:
            continue

        node_cutoff_score[node] = edges_to_keep[-1][2]["Cosine"]

        #print("DELETE", edges_to_delete)


        #for edge_to_remove in edges_to_delete:
        #    G.remove_edge(edge_to_remove[0], edge_to_remove[1])


    #print("After Top K", len(G.edges()))
    #Doing this for each pair, makes sure they are in each other's top K
    edge_to_remove = []
    for edge in G.edges(data=True):
        cosine_score = edge[2]["Cosine"]
        threshold1 = node_cutoff_score[edge[0]]
        threshold2 = node_cutoff_score[edge[1]]

        if cosine_score < threshold1 or cosine_score < threshold2:
            edge_to_remove.append(edge)

    for edge in edge_to_remove:
        G.remove_edge(edge[0], edge[1])

    #print("After Top K Mutual", len(G.edges()))