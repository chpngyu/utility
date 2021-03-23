import sys, os
# for increasing maximum recursion of heapmap
#sys.setrecursionlimit(100000)

import os.path as op
def basename(path):### get base name of a path by removing its extension
    name = op.basename(path)
    if '.' in name:
        name = name.split('.')[0]
    return name

import pandas as pd
def unit_expression(data_frame):
    columns = data_frame.columns
    for index in data_frame.index:
        values = data_frame.ix[index].values
        max_v, min_v = max(values), min(values)
        data_frame.loc[index] = [(v-min_v)/(max_v - min_v) for v in values]

from scipy import stats
def zscore_expression(data_frame):
    for index in data_frame.index:
        data_frame.ix[index] = stats.zscore(data_frame.ix[index])

def my_correlation(x, y):
    n = len(x)
    sum_x = float(sum(x))
    sum_y = float(sum(y))
    sum_x_sq = sum(map(lambda x: pow(x, 2), x))
    sum_y_sq = sum(map(lambda x: pow(x, 2), y))
    psum = sum(map(lambda x, y: x * y, x, y))
    num = psum - (sum_x * sum_y/n)
    den = pow((sum_x_sq - pow(sum_x, 2) / n) * (sum_y_sq - pow(sum_y, 2) / n), 0.5)
    if den == 0:
        return 1
    dist = 1 - num/den
    return max(dist, 0)

import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from sklearn.cluster import AgglomerativeClustering
from random import shuffle

parser = argparse.ArgumentParser(description='Plot Heatmap')
parser.add_argument('-c', '--cluster', type=str, required=True,
                    help = 'Number of clusters or provide o file that contains category for each item')
parser.add_argument('-d', '--data', type=str, required=True,
                    help='A data matrix')
parser.add_argument('-f', '--format', type=int, default=0, choices=[0,1],
                    help='data format in CSV (0: default) or Tab (1)')
parser.add_argument('-i', '--index', type=str, required=True,
                    help='Index column header in the given data matrix. The value in the index column must be string not number')
parser.add_argument('--filter', type=str,
                    help='Only use the given list of elements to generate heatmap')
parser.add_argument('-m', '--metric', type=str, default='correlation',
                    help = 'Distance function (default: correlation)')
parser.add_argument('-x', '--xcluster', action='store_true',
                    help = 'cluster columns (default: False)')
parser.add_argument('-y', '--ycluster', action='store_true',
                    help = 'cluster row (default: True)')
parser.add_argument('--labelx', action='store_true',
                    help = 'Show x-tick labels')
parser.add_argument('--labely', action='store_true',
                    help = 'Show y-tick labels')
parser.add_argument('--minstd', type=float, default=0,
                    help = 'Retain elements which have std > minstd (default: 0)')
parser.add_argument('--scale', type=int, default=0,
                    help = 'Scaling expression levels to unvaried(0), unit(1) or z-score(2). Default is unvaried(0)')
parser.add_argument('--mean', type=str,
                    help = 'Output mean profile for each cluster (option)')

args = parser.parse_args()

if args.format == 0:
    x = pd.read_csv(args.data, index_col=args.index, sep=',')
else:
    x = pd.read_csv(args.data, index_col=args.index, sep='\t')

if args.filter != None:
    deg = open(args.filter, 'U').read().split()
    print('Got %d elements from %s' % (len(deg), args.filter))
    deg = set(deg) & set(x.index)
    print('where %d elements are common with %s' % (len(deg), args.data))
    x = x.loc[deg]
    #x.to_csv("~/Downloads/saved.peak.counts.txt",sep='\t')
    print('Got %d elements after filtering' % len(x))
    if len(x) < 2:
        print('Only one or no element remained. Program terminated!')
        exit(0)

x = x.loc[x.std(axis=1)>args.minstd]
print('%d elements are retained with std > %g' % (len(x), args.minstd))
if len(x) < 2:
    print('Only one gene remained. Program terminated!')
    exit(0)
if args.scale:
    scale_expression = [unit_expression, zscore_expression][args.scale-1]
    scale_func = ['unit', 'z-score'][args.scale-1]
    print('Scaling expression set to %s' % scale_func)
    scale_expression(x)


# assign color for each item (row) by using
if os.path.isfile(args.cluster): # user provided clusters
   item2categ = {}
   with open(args.cluster, 'U') as ifile:
       for line in ifile:
           line = line.strip()
           if len(line) == 0:
               continue
           item, categ = line.split('\t')
           item2categ[item] = categ
   category = set(item2categ.values())
   no_cluster = len(category)
   print('Got %d items with totally %d categories from %s' % (len(item2categ), no_cluster, args.cluster))
   cluster_pal = sns.hls_palette(no_cluster, l=.4, s=1)
   shuffle(cluster_pal)
   category2index = dict(zip(category, range(no_cluster)))
   y_pred = [category2index[item2categ[item]] for item in x.index]
   cluster_colors = [cluster_pal[y] for y in y_pred]
else: # agglomerative clustering
    no_cluster = int(args.cluster)
    y_pred = AgglomerativeClustering(n_clusters=no_cluster, affinity=args.metric, linkage='average').fit_predict(x.values)
    cluster_pal = sns.hls_palette(no_cluster, l=.4, s=1)
    shuffle(cluster_pal)
    cluster_pal = dict(zip(range(no_cluster), cluster_pal))
    cluster_colors = [cluster_pal[idx] for idx in y_pred]

if args.metric == 'correlation':
    metric = my_correlation
else:
    metric = args.metric
cmap = sns.diverging_palette(h_neg=250, h_pos=10, s=99, l=50, as_cmap=True)
#cmap = sns.palplot(sns.cubehelix_palette(10, start=2, rot=0, dark=0, light=.95, reverse=True))
#cmap = sns.diverging_palette(h_neg=130, h_pos=10, s=99, l=55, center='dark', as_cmap=True)
cg = sns.clustermap(x, metric=metric, cmap=cmap, #row_colors=cluster_colors,
                    xticklabels=args.labelx, yticklabels=args.labely,
                    col_cluster=args.xcluster, row_cluster=args.ycluster, figsize=(6, 14))

cg.ax_heatmap.set_xticklabels(cg.ax_heatmap.get_xmajorticklabels(), fontsize = 9)
cg.ax_heatmap.set_yticklabels(cg.ax_heatmap.get_ymajorticklabels(), fontsize = 9)
plt.setp(cg.ax_heatmap.xaxis.get_majorticklabels(), rotation=90)
plt.setp(cg.ax_heatmap.yaxis.get_majorticklabels(), rotation=0)

filename = 'Heatmap_' + basename(args.data) + '.pdf'
plt.show()
#plt.savefig(filename, format='pdf')
print('Saving figure to %s' % filename)
filename = 'Clusters_' + basename(args.data)
ofile = open(filename, 'w')
y_cluster_order = []
if args.ycluster:
    ylabels = [a._text for a in cg.ax_heatmap.get_yticklabels()]
    y_pred_dict = {str(u):v for u, v in zip(x.index, y_pred)}
    for idx in range(len(ylabels)-1, -1, -1):
        y = str(ylabels[idx])
        ofile.write('%s\t%d\n' % (y, y_pred_dict[y]))
        if y_pred_dict[y] not in y_cluster_order:
            y_cluster_order.append(y_pred_dict[y])
else:
    for idx, row in enumerate(x.iterrows()):
        ofile.write('%s\t%d\n' % (row[0], y_pred[idx]))
        if y_pred[idx] not in y_cluster_order:
            y_cluster_order.append(y_pred[idx])
ofile.close()
print('Ouputing cluster info to %s' % filename)
if args.xcluster & args.labelx:
    xticks_in_order = [xlabels._text for xlabels in cg.ax_heatmap.get_xticklabels()]
    x = x[xticks_in_order]
    xticks_in_order = '\n'.join(xticks_in_order)
    print('Clustering sample:\n% s' % xticks_in_order)

if args.mean != None:
    def get_mean(data_frame):
        _mean_profile = []
        for col in data_frame.columns:
            _mean_profile.append(data_frame[col].mean())
        return _mean_profile
    cluster_data = []
    inedx_col = []
    cluster_size = []
    for c in y_cluster_order:
        genes = [u for u in ylabels if y_pred_dict[u] == c]
        mprofile = get_mean(x.loc[genes])
        cluster_data.append(mprofile)
        inedx_col.append('C%d' % c)
        cluster_size.append(len(genes))
    cluster_data = pd.DataFrame(cluster_data, index=inedx_col, columns=x.columns)
    plt.clf()
    heatmap = sns.heatmap(cluster_data, cmap=cmap)
    plt.setp(heatmap.xaxis.get_majorticklabels(), rotation=90)
    plt.savefig(args.mean)
    print('A figure of mean profiles has been written to %s' % args.mean)
    cluster_data['Size'] = cluster_size
    cluster_data.to_csv('mean_cluster_size.txt', sep='\t')
