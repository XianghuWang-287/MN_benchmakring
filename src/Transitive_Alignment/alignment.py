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
import math
import time
import os
import argparse

SpectrumTuple = collections.namedtuple(
    "SpectrumTuple", ["precursor_mz", "precursor_charge", "mz", "intensity"]
)
def _cosine_fast(
    spec: SpectrumTuple,
    spec_other: SpectrumTuple,
    fragment_mz_tolerance: float,
    allow_shift: bool,
) -> Tuple[float, List[Tuple[int, int]]]:
    precursor_charge = max(spec.precursor_charge, 1)
    precursor_mass_diff = (
        spec.precursor_mz - spec_other.precursor_mz
    ) * precursor_charge
    # Only take peak shifts into account if the mass difference is relevant.
    num_shifts = 1
    if allow_shift and abs(precursor_mass_diff) >= fragment_mz_tolerance:
        num_shifts += precursor_charge
    other_peak_index = np.zeros(num_shifts, np.uint16)
    mass_diff = np.zeros(num_shifts, np.float32)
    for charge in range(1, num_shifts):
        mass_diff[charge] = precursor_mass_diff / charge

    # Find the matching peaks between both spectra.
    peak_match_scores, peak_match_idx = [], []
    for peak_index, (peak_mz, peak_intensity) in enumerate(
        zip(spec.mz, spec.intensity)
    ):
        # Advance while there is an excessive mass difference.
        for cpi in range(num_shifts):
            while other_peak_index[cpi] < len(spec_other.mz) - 1 and (
                peak_mz - fragment_mz_tolerance
                > spec_other.mz[other_peak_index[cpi]] + mass_diff[cpi]
            ):
                other_peak_index[cpi] += 1
        # Match the peaks within the fragment mass window if possible.
        for cpi in range(num_shifts):
            index = 0
            other_peak_i = other_peak_index[cpi] + index
            while (
                other_peak_i < len(spec_other.mz)
                and abs(
                    peak_mz - (spec_other.mz[other_peak_i] + mass_diff[cpi])
                )
                <= fragment_mz_tolerance
            ):
                peak_match_scores.append(
                    peak_intensity * spec_other.intensity[other_peak_i]
                )
                peak_match_idx.append((peak_index, other_peak_i))
                index += 1
                other_peak_i = other_peak_index[cpi] + index

    score, peak_matches = 0.0, []
    if len(peak_match_scores) > 0:
        # Use the most prominent peak matches to compute the score (sort in
        # descending order).
        peak_match_scores_arr = np.asarray(peak_match_scores)
        peak_match_order = np.argsort(peak_match_scores_arr)[::-1]
        peak_match_scores_arr = peak_match_scores_arr[peak_match_order]
        peak_match_idx_arr = np.asarray(peak_match_idx)[peak_match_order]
        peaks_used, other_peaks_used = set(), set()
        for peak_match_score, peak_i, other_peak_i in zip(
            peak_match_scores_arr,
            peak_match_idx_arr[:, 0],
            peak_match_idx_arr[:, 1],
        ):
            if (
                peak_i not in peaks_used
                and other_peak_i not in other_peaks_used
            ):
                score += peak_match_score
                # Save the matched peaks.
                peak_matches.append((peak_i, other_peak_i))
                # Make sure these peaks are not used anymore.
                peaks_used.add(peak_i)
                other_peaks_used.add(other_peak_i)

    return score, peak_matches

def peak_tuple_to_dic(peakmatches):
    dic={}
    for peakmatch in peakmatches:
        dic[peakmatch[0]]=peakmatch[1]
    return dic

def norm_intensity(intensity):
    return np.copy(intensity)/np.linalg.norm(intensity)
    #return(intensity)

def realign_path(path):
    final_match_list=[]
    spec_1 = spec_dic[path[0]]
    spec_2 = spec_dic[path[0+1]]
    score,peak_matches = _cosine_fast(spec_1,spec_2,0.5,True)
    final_match_list=peak_matches
    idx=1
    while (len(final_match_list)!=0 and idx <len(path)-1):
        temp_peakmatch_list=[]
        spec_1 = spec_dic[path[idx]]
        spec_2 = spec_dic[path[idx+1]]
        score,peak_matches = _cosine_fast(spec_1,spec_2,0.5,True)
        peak_dic1=peak_tuple_to_dic(final_match_list)
        peak_dic2=peak_tuple_to_dic(peak_matches)
        for key, value in peak_dic1.items():
            if (peak_dic2.get(value)):
                temp_peakmatch_list.append((key,peak_dic2[value]))
        final_match_list=temp_peakmatch_list
        idx=idx+1
    spec_start = spec_dic[path[0]]
    spec_end = spec_dic[path[-1]]
    _, matched_peaks = _cosine_fast(spec_start,spec_end,0.5,True)
    peak_match_scores = []
    intesity1=spec_dic[path[0]][3]
    intesity2=spec_dic[path[-1]][3]
    if (len(final_match_list)):
        final_match_list=final_match_list+matched_peaks
        for matched_peak in final_match_list:
            peak_match_scores.append(intesity1[matched_peak[0]]*intesity2[matched_peak[1]])
    score, peak_matches = 0.0, []
    if len(peak_match_scores) > 0:
        # Use the most prominent peak matches to compute the score (sort in
        # descending order).
        peak_match_scores_arr = np.asarray(peak_match_scores)
        peak_match_order = np.argsort(peak_match_scores_arr)[::-1]
        peak_match_scores_arr = peak_match_scores_arr[peak_match_order]
        peak_match_idx_arr = np.asarray(final_match_list)[peak_match_order]
        peaks_used, other_peaks_used = set(), set()
        for peak_match_score, peak_i, other_peak_i in zip(
            peak_match_scores_arr,
            peak_match_idx_arr[:, 0],
            peak_match_idx_arr[:, 1],
        ):
            if (
                peak_i not in peaks_used
                and other_peak_i not in other_peaks_used
            ):
                score += peak_match_score
                # Save the matched peaks.
                peak_matches.append((peak_i, other_peak_i))
                # Make sure these peaks are not used anymore.
                peaks_used.add(peak_i)
                other_peaks_used.add(other_peak_i)
        return peak_matches, score
    else:
        return "no match",_
def re_alignment_parallel(args):
    node1,node2=args
    if nx.has_path(G_all_pairs, node1, node2):
        key_path_start_time = time.time()
        all_shortest_hops = [p for p in nx.all_shortest_paths(G_all_pairs, node1, node2, weight=None, method='dijkstra')]
        Max_path_weight = 0
        Max_path = []
        for hop in all_shortest_hops:
            Path_weight = 0
            for node_i in range(len(hop) - 1):
                Path_weight = Path_weight + G_all_pairs[hop[node_i]][hop[node_i + 1]]['Cosine']
            if (Path_weight > Max_path_weight):
                Max_path_weight = Path_weight
                Max_path = hop
        key_path_end_time = time.time()
        key_path_running_time = key_path_end_time - key_path_start_time
        align_peak_start_time = time.time()
        matched_peaks,score = realign_path(Max_path)
        align_peak_end_time = time.time()
        align_peak_running_time = align_peak_end_time - align_peak_start_time
        if score >= 0.85:
            print(Max_path)
        if matched_peaks != "no match":
            G_all_pairs_realignment.add_edge(node1,node2)
            G_all_pairs_realignment[node1][node2]['Cosine']=score
            return (node1,node2,score,key_path_running_time,align_peak_running_time)
        else:
            return
    else:
        return


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
        print("starting align library:"+library)
        summary_file_path = "../../data/summary/"+library+"_summary.tsv"
        merged_pairs_file_path = "../../data/merged_pairs/"+library+"_merged_pairs.tsv"
        mgf_file_path = "../../data/converted/"+library+"_converted.mgf"
        cluster_summary_df = pd.read_csv(summary_file_path)
        all_pairs_df = pd.read_csv(merged_pairs_file_path, sep='\t')
        G_all_pairs = nx.from_pandas_edgelist(all_pairs_df, "CLUSTERID1", "CLUSTERID2", "Cosine")
        G_all_pairs_realignment = G_all_pairs.copy()
        spec_dic = {}
        print("Start create spectrum dictionary")
        for spectrum in tqdm(mgf.read(mgf_file_path)):
            params = spectrum.get('params')
            precursor_mz = cluster_summary_df.loc[int(params['scans']) - 1]["Precursor_MZ"]
            charge = cluster_summary_df.loc[int(params['scans']) - 1]["Charge"]
            mz_array = spectrum.get('m/z array')
            intensity_array = spectrum.get('intensity array')
            filtered_mz = []
            filtered_intensities = []
            precursor_value = float(cluster_summary_df.loc[cluster_summary_df['scan'] == int(params['scans'])]["Precursor_MZ"].values[0])
            for i, mz in enumerate(mz_array):
                peak_range = [j for j in range(len(mz_array)) if abs(mz_array[j] - mz) <= 25]
                sorted_range = sorted(peak_range, key=lambda j: intensity_array[j], reverse=True)
                if i in sorted_range[:6]:
                    if abs(mz - precursor_value) > 17:
                        filtered_mz.append(mz)
                        filtered_intensities.append(intensity_array[i])
            filtered_intensities = [math.sqrt(x) for x in filtered_intensities]
            spec_dic[int(params['scans'])] = SpectrumTuple(precursor_value, charge, filtered_mz, norm_intensity(filtered_intensities))

        with Pool(processes=28, maxtasksperchild=1000) as pool:

            values = [[node1, node2] for [node1, node2] in nx.non_edges(G_all_pairs)]

            results = list(tqdm(pool.imap(re_alignment_parallel, values), total=len(values)))
            # print the results
        key_path_running_list = [result[3] for result in results]
        align_peak_running_list = [result[4] for result in results]
        print("key path total running time:", sum(key_path_running_list))
        print("align peak running time:", sum(align_peak_running_list))
        result_file_path = "../../results/alignment_results/" + library + "_realignment.pkl"
        with open(result_file_path, 'wb') as f:
            pickle.dump(results, f)
    print("Finish align all libraries")

