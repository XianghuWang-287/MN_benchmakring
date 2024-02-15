from Classic import *
import pickle
import argparse
import random

def add_edges_to_mst(original_graph, mst):
    remaining_edges = [(u, v, original_graph[u][v]['Cosine']) for u, v in original_graph.edges() if not mst.has_edge(u, v)]
    remaining_edges.sort(key=lambda x: x[2], reverse=True)

    average_weight = calculate_average_weight(mst)

    for u, v, weight in remaining_edges:
        mst.add_edge(u, v, Cosine=weight)
        new_average_weight = calculate_average_weight(mst)
        if new_average_weight >= average_weight:
            mst.remove_edge(u, v)
            break
        average_weight = new_average_weight

    return mst

def calculate_average_weight(graph):
    total_weight = sum(graph[u][v]['Cosine'] for u, v in graph.edges())
    average_weight = total_weight / graph.number_of_edges()
    return average_weight

def polish_subgraph(G):
    if G.number_of_edges() == 0:
        return G
    maximum_spanning_tree = nx.maximum_spanning_tree(G, weight='Cosine')
    polished_subgraph = add_edges_to_mst(G, maximum_spanning_tree)
    return polished_subgraph

def sample_graph_by_probability(graph, sample_percentage):
    node_list=[]
    for node in graph.nodes():
        if random.random() <= sample_percentage:
            node_list.append(node)
    return graph.subgraph(node_list).copy()


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
        summary_file_path = "../../data/summary/"+library+"_summary.tsv"
        merged_pairs_file_path = "../../data/merged_pairs/"+library+"_merged_pairs.tsv"
        re_align_edge_file_path = "../../results/alignment_results/"+library+"_realignment.pkl"
        cluster_summary_df = pd.read_csv(summary_file_path)
        all_pairs_df = pd.read_csv(merged_pairs_file_path, sep='\t')
        G_all_pairs = nx.from_pandas_edgelist(all_pairs_df, "CLUSTERID1", "CLUSTERID2", "Cosine")
        print('graph with {} nodes and {} edges'.format(G_all_pairs.number_of_nodes(), G_all_pairs.number_of_edges()))
        print("constructing dic for finger print")
        dic_fp = fingerprint_dic_construct_InCHI(cluster_summary_df)
        print("constructing the re-alignment graph")
        G_all_pairs_realignment = G_all_pairs.copy()
        with open(re_align_edge_file_path, 'rb') as f:
            realignment_edgelist = pickle.load(f)
        new_list=[]
        for item in realignment_edgelist:
            if(item != None):
                if item[2]>=0.75:
                    new_list.append(item)
        for item in new_list:
            if (item != None):
                G_all_pairs_realignment.add_edge(item[0], item[1], Cosine=item[2])
        results_df_list = []
        thresholds = [x / 100 for x in range(70, 95)]
        for threshold in tqdm(thresholds):
            cast_cluster = CAST_cluster(G_all_pairs_realignment, threshold)
            cast_score_list = []
            cast_components = [G_all_pairs_realignment.subgraph(c).copy() for c in cast_cluster]
            benchmark_set = []
            for component in cast_components:
                benchmark_set.append(polish_subgraph(component))
            for component in benchmark_set:
                cast_score_list.append(subgraph_score_dic(component, cluster_summary_df, dic_fp))
            cast_number = [len(x) for x in cast_cluster]
            df_cast = pd.DataFrame(list(zip(cast_score_list, cast_number)), columns=['score', 'number'])
            results_df_list.append(df_cast)
        result_file_path = "../../results/results-re-cast/" + library + "_re_cast_benchmark.pkl"
        with open(result_file_path, 'wb') as file:
            pickle.dump(results_df_list, file)





