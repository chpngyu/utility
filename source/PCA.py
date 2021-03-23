####read a data table and perform PCA analysis
#
# global parameters
delta_shift = 0.01
sample_label_size = 4
tick_label_size = 10
axis_label_size = 14
titile_font_size = 16
point_size = 50
alphav = 0.7
legend_fontsize=8
legend_location=2
tick_pad = 2 # space between the ticklabels and the axes

def unit_expression(data_frame):
    for index in data_frame.index:
        values = data_frame.ix[index].values
        max_v, min_v = max(values), min(values)
        data_frame.loc[index] = [(v-min_v)/(max_v - min_v) for v in values]

import matplotlib
from matplotlib import pyplot as plt
from sklearn import decomposition # for PCA
import numpy as np
import pandas as pd

from mpl_toolkits.mplot3d import Axes3D

def find_category(sample_names, category_file):
    record_group = {}
    with open(category_file) as ifile:
        for line in ifile:
            line = line.strip()
            if len(line) == 0:
                continue
            record, group = line.split('\t')
            record_group[record] = group
    _groups = sorted(set(record_group.values()))
    print('Got %d records with %d categories from %s' % (len(record_group), len(_groups), category_file))
    _sample_groups = []
    for sample in sample_names:
        idx = _groups.index(record_group[str(sample)])
        _sample_groups.append(idx)
    return _sample_groups, dict([(_idx, u) for _idx, u in enumerate(_groups)])

def data_in_category(_data, _sample_groups, _category):
    for catg in _category.keys():
        _result = ()
        for idx, _group in enumerate(_sample_groups):
            if _group == catg:
                _result += (_data[idx],)
        if len(_result):
            yield np.vstack(_result), catg

import argparse
parser = argparse.ArgumentParser(description='PCA analysis')
parser.add_argument('-d', '--data', type=str, required=True, help = 'input data')
parser.add_argument('-f', '--format', type=int, default=0, choices=[0,1],
                    help='data format in CSV (0: default) or Tab (1)')
parser.add_argument('-b', '--box', type=str,
                    help='output file for box plot')
parser.add_argument('-i', '--index', type=str, help='provide index column (default: first column)')
parser.add_argument('--dimension', type=int, default=2, choices=[2,3],
                    help='how many dimentions to plot the result of PCA (default: 2)')
parser.add_argument('--way', type=int, choices=[0,1], default=1,
                    help='use row (0) or column (1) to compute PCA (default: 1 (column))')
parser.add_argument('--filter', type=str,
                    help='Only use the given list of records in row to generate heatmap')
parser.add_argument('--scale', action='store_true',
                    help = 'Scaling expression levels to be unit')
parser.add_argument('--log', action='store_true',
                    help='taking log(x+1) for expression')
parser.add_argument('--save', type=str,
                    help = 'save figure to a file (if not set, figure will display by GUI)')
parser.add_argument('--pcav', type=str,
                    help = 'write PCA values to a file')
parser.add_argument('-c', '--category', type=str,
                    help = 'category for each sample (option)')
parser.add_argument('--label', nargs='?', const='all',
                    help = 'show label for each sample (default=No)')



args = parser.parse_args()
if args.index == None:
    index_col = 0
else:
    index_col = args.index
if args.format:
    data = pd.read_csv(args.data, index_col=index_col, sep='\t')
else:
    data = pd.read_csv(args.data, index_col=index_col)

if args.filter != None:
    deg = open(args.filter).read().split()
    data = data.loc[deg]
    print('Got %d genes after filtering' % len(data))

print('Size: %d*%d' % (len(data), len(data.columns)))
data.dropna(inplace=True)
print('After dropping NA, size becomes: %d*%d' % (len(data), len(data.columns)))

if args.way:
    print('Using columns as samples to perfrom PCA')
    data = data.transpose()
sampleIDs = data.index

if args.scale:
    print('Scaling values to be unit')
    unit_expression(data)

if args.log:
    import numpy
    print('Taking log for values')
    data = data+1
    data = data.apply(numpy.log)

matplotlib.rc('xtick', labelsize=tick_label_size)
matplotlib.rc('ytick', labelsize=tick_label_size)
plt.rcParams['xtick.major.pad']= ('%d' % tick_pad)
plt.rcParams['ytick.major.pad']= ('%d' % tick_pad)
if args.box:
    data.boxplot(rot=18, fontsize=12, return_type='axes')
    plt.xticks(rotation='vertical')
    plt.subplots_adjust(bottom=0.15)
    plt.savefig(args.box)
    print('Box plot has been saved to %s' % args.box)
    plt.clf()

pca = decomposition.PCA()
pca.fit(data.values)
results = pca.transform(data)
print('Component\tExplained variance')
for idx, var in enumerate(pca.explained_variance_ratio_):
    if idx > 10:
        print('...')
        break
    else:
        print('%d\t%g' % (idx+1, var))

if args.category != None:
    sample_groups, groups = find_category(data.index, args.category)
    no_group = len(groups)
    if no_group <= 10:
        colors = plt.cm.tab10(np.linspace(0, 1, no_group))
    else:
        colors = plt.cm.hsv(np.linspace(0, 1, no_group))
    print('Category\tNo. of items')
    if args.dimension == 2:
        for grp_data, label in data_in_category(results, sample_groups,  groups):
            color = list(colors[label])
            color[-1] = alphav
            no_data = len(grp_data)
            color = np.stack(color*no_data).reshape(no_data,4)
            print('%s\t%d' % (groups[label], no_data))
            plt.scatter(grp_data[:,0], grp_data[:,1], s=point_size,
                        c=color, label=groups[label])
    else:
        plt.rcParams['axes.facecolor']='white'
        fig = plt.figure(1)
        ax = Axes3D(fig, rect=[0, 0, .95, 1])
        for grp_data, label in data_in_category(results, sample_groups,  groups):
            color = list(colors[label])
            color[-1] = alphav
            no_data = len(grp_data)
            color = np.stack(color*no_data).reshape(no_data,4)
            print('%s\t%d' % (groups[label], len(grp_data)))
            ax.scatter(grp_data[:,0], grp_data[:,1], grp_data[:,2], s=point_size,
                       c=color, label=groups[label])
    plt.legend(fontsize=legend_fontsize, loc=legend_location)
else:
    if args.dimension == 2:
        plt.scatter(results[:,0], results[:,1], s=point_size, alpha=alphav)
    else:
        plt.rcParams['axes.facecolor']='white'
        fig = plt.figure(1)
        ax = Axes3D(fig, rect=[0, 0, .95, 1])
        ax.scatter(results[:,0], results[:,1], results[:,2], s=point_size)
plt.title('Principal component analysis', fontsize=titile_font_size)
label_value = 'PC1 ({:.1f}%)'.format(pca.explained_variance_ratio_[0]*100)
plt.xlabel(label_value, fontsize=axis_label_size)
label_value = 'PC2 ({:.1f}%)'.format(pca.explained_variance_ratio_[1]*100)
plt.ylabel(label_value, fontsize=axis_label_size)

if args.dimension == 3:
    label_value = 'PC3 ({:.1f}%)'.format(pca.explained_variance_ratio_[2]*100)
    ax.set_zlabel(label_value, fontsize=axis_label_size)
# labeling sample IDs
if args.label != None:
    delta_xpos = (max(results[:,0])-min(results[:,0]))*delta_shift
    delta_ypos = (max(results[:,1])-min(results[:,1]))*delta_shift
    if args.label != 'all':
        used_records = [x.strip() for x in open(args.label, 'U')]
        print('Got %d records from %s and will for labeling' % (len(used_records), args.label))
    if args.dimension == 2:
        for ix, pos in enumerate(zip(results[:,0], results[:,1])):
            if (args.label != 'all') and (sampleIDs[ix] not in used_records):
                continue
            plt.text(pos[0]+delta_xpos, pos[1]+delta_ypos, sampleIDs[ix], fontsize=sample_label_size)
    else:
        delta_zpos = (max(results[:,2])-min(results[:,2]))*delta_shift
        for ix, pos in enumerate(zip(results[:,0], results[:,1], results[:,2])):
            if (args.label != 'all') and (sampleIDs[ix] not in used_records):
                continue
            ax.text(pos[0]+delta_xpos, pos[1]+delta_ypos, pos[2]+delta_zpos, sampleIDs[ix], fontsize=sample_label_size)

if args.save != None:
    if args.dimension == 2:
        plt.savefig(args.save, bbox_inches='tight', pad_inches=0.3)
    else:
        plt.savefig(args.save, facecolor=fig.get_facecolor(), bbox_inches='tight', pad_inches=0.3)
    print('Figure has been saved to %s' % args.save)
else:
    plt.show()

if args.pcav != None:
    pca_values = {}
    pca_values['PCA1'] = results[:,0]
    pca_values['PCA2'] = results[:,1]
    if args.dimension == 3:
        pca_values['PCA3'] = results[:,2]
    pd.DataFrame(pca_values, index=data.index).to_csv(args.pcav)
    print('PCA values for each sample have been written to ' + args.pcav)