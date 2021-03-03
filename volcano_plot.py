# -*- coding: utf-8 -*-
"""
Created on Mar 6, 2021

@author: Chun-Ping Yu
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
parser = argparse.ArgumentParser('Volcano plot')
# add options
parser.add_argument('-d', '--data', type=str, required=True,
					help='input data')
parser.add_argument('--sep', type=str, default=',',
					help='Delimiter-separated values for the given data (default: ,)')
parser.add_argument('--index', type=int, default=0,
					help='index column (default: 0)')
parser.add_argument('--fc', type=int, default=1,
					help='the number of column for log2 fold change (default: 1)')
parser.add_argument('--pvalue', type=int, default=2,
					help='the number of column for pvalue (default: 2)')
parser.add_argument('--threshold', type=float, default=0.05,
					help='highlight the points that have pvalue less thant it (default: 0.05)')
parser.add_argument('--highlight', type=str, default='red',
					help='set color for the highlight points (default: red)')
parser.add_argument('--color', type=str, default='gray',
					help='set color for the other points with p-value >= threshold (default: gray)')
parser.add_argument('--title', type=str,
					help='set title of the plot if give')
parser.add_argument('-o', '--output', type=str, default='volcano_plot.jpg',
					help='filename for the plot (volcano_plot.jpg)')
# parse arguments
args = parser.parse_args()

data = pd.read_csv(args.data, index_col=args.index)
print(f'Read data from {args.output} with size {data.shape}')
fc = data.columns[args.fc-1]
pvalue = data.columns[args.pvalue-1]

filter = data[pvalue] < args.threshold
background = data[~filter]
highlight = data[filter]

plt.scatter(background[fc].values, -np.log10(background[pvalue].values), color=args.color, alpha=0.6)
plt.scatter(highlight[fc].values, -np.log10(highlight[pvalue].values), color=args.highlight)
plt.xlabel(fc)
plt.ylabel(f'-log10 {pvalue}')
if args.title:
	plt.title(args.title)
plt.tight_layout()
plt.savefig(args.output)
print(f'Done. Please see {args.output}')