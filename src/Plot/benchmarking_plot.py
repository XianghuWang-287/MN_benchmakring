
import matplotlib.pyplot as plt
import numpy as np
import pickle


def weighted_average(dataframe, value, weight):
    val = dataframe[value]
    wt = dataframe[weight]
    return (val * wt).sum() / wt.sum()

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

with open("../../results/results-classic/GNPS-NIH-SMALLMOLECULEPHARMACOLOGICALLYACTIVE_classic_benchmark.pkl","rb") as f:
    results_df_list = pickle.load(f)
with open("../../results/results-re-cast/GNPS-NIH-SMALLMOLECULEPHARMACOLOGICALLYACTIVE_re_cast_benchmark.pkl","rb") as f:
    results_re_MST_df_list = pickle.load(f)
with open("../../results/results-cast/GNPS-NIH-SMALLMOLECULEPHARMACOLOGICALLYACTIVE_cast_benchmark.pkl","rb") as f:
    results_cast_df_list = pickle.load(f)

x_default = 7
y_default = 0.69670551
X_results_base=[1175,511,3,1,1,1]
y_results_base=[0.4254776,0.49974456,0.64723446,0.77786232,0.86408052,0.89820489]
x_stru_results = [1186,1112,804,309,19,10,6,5,3,3,2]
y_stru_results = [0.50579754,0.59545469,0.72411139,0.80563852,0.86268481,0.90332193,0.93423447,0.95666362,0.97475872,0.98705394,0.99585878]
x_ms2ds_fixed_results_base = [1460.0, 1460.0, 1460.0, 1460.0, 1429.0, 1408.0, 1364.0, 1285.0, 1123.0, 831.0, 521.0, 8.0, 2.0, 1, 1, 1, 1, 1, 1]
y_ms2ds_fixed_results_base = [0.36341719449667564, 0.3627775598384254, 0.36235338108661874, 0.3651573467458511, 0.40121577864573016, 0.4131286395487256, 0.43262271978716205, 0.4664619474714032, 0.5146205082208899, 0.5948417150119151, 0.6746627633268013, 0.7344754177836268, 0.8416950909957551, 0.9149834914386169, 0.9137305005029621, 0.9146864703726121, 0.9184357158445278, 0.9206985762347929, 0.954136960740777]

fig = plt.figure(figsize=(12,8))
x_results = np.array([cal_N50(x,1460, 0.25) for x in results_df_list])
y_results = np.array([weighted_average(x, 'score', 'number') for x in results_df_list])
x_re_MST_results = np.array([cal_N50(x,1201, 0.18) for x in results_re_MST_df_list])
y_re_MST_results = np.array([weighted_average(x, 'score', 'number') for x in results_re_MST_df_list])
x_cast_results = np.array([cal_N50(x,1201, 0.3) for x in results_cast_df_list])
y_cast_results = np.array([weighted_average(x, 'score', 'number') for x in results_cast_df_list])



plt.scatter(X_results_base, y_results_base,s=50,color='brown', label = 'Baseline')
plt.scatter(x_results, y_results,s=50, label = 'Classic')
plt.scatter(x_default, y_default,s=200,marker = '*', color = 'orange',label = 'Classic-default')
plt.scatter(x_ms2ds_fixed_results_base, y_ms2ds_fixed_results_base,s=50,color = 'purple', label = 'MS2DeepScore')
plt.scatter(x_cast_results, y_cast_results,marker = '^',s=50,color = 'limegreen', label = 'CAST')
plt.scatter(x_re_MST_results, y_re_MST_results,s=50, marker = ',',color='r',label = 'CAST + Transitive Alignment')
plt.scatter(x_stru_results,y_stru_results,s=50,label = 'Optimal Network Topology')

plt.xlabel('N20',fontsize = 20)
plt.ylabel('Network Accuracy Score',fontsize = 20)
plt.xticks(fontsize = 20)
plt.yticks(fontsize = 20)
plt.xscale('log')
plt.legend(loc='lower left',prop={'size': 14})
plt.title('Benchmarking results for NIH SPAC', fontsize= 22)
plt.show()