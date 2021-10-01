##this script confirms our implementation of the Manczinger et al. algorithm reproduces their values

import pandas as pd
import pylab as plt
import numpy as np
from kl import kl

import pathlib
path = pathlib.Path.cwd()
if path.stem == 'iedb':
    cwd = path
else:
    cwd = list(path.parents)[::-1][path.parts.index('iedb')]
    import sys
    sys.path.append(str(cwd))


##get original kl values
fig1_data = pd.read_csv(cwd / 'data' / 'fig1_source_data.csv') #from https://www.nature.com/articles/s43018-021-00226-4
original = {'HLA-' + i[0] + '*' + i[1:]: j for i, j in zip(fig1_data['Allele name'].values, fig1_data['Allele promiscuity'].values)}

##the peptides they used
mhc_df = pd.read_csv(cwd / 'data' / 'mhc_ligand.txt', sep='\t') #from https://www.nature.com/articles/s43018-021-00226-4
mhc_df['Description'] = mhc_df['Description'].str.upper()
mhc_df['Allele Name'] = mhc_df['Allele Name'].str.replace(':', '')

##restrict data to peptide lengths used
iedb_df = mhc_df.loc[(mhc_df['Description'].str.len() < 13) & (mhc_df['Description'].str.len() > 7)]

prom_dict = {}
for allele in mhc_df['Allele Name'].unique():
    if sum(mhc_df['Allele Name'] == allele) >= 400:
        prom_dict[allele] = kl(mhc_df.loc[mhc_df['Allele Name'] == allele])[0]


from matplotlib import cm
paired = [cm.get_cmap('Paired')(i) for i in range(12) if i not in [4, 5]]
fig = plt.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(top=0.985,
bottom=0.1,
left=0.09,
right=0.975,
hspace=0.2,
wspace=0.2)

ax.scatter([original[i] for i in original], [prom_dict[i] for i in original], color=paired[1], edgecolor='k')
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_linewidth(1)
ax.spines['left'].set_position(['outward', 2])
ax.spines['left'].set_bounds(1.4, 3.5)
ax.spines['bottom'].set_linewidth(1)
ax.spines['bottom'].set_bounds(1.4, 3.5)
ax.spines['bottom'].set_position(['outward', 2])
ax.set_xticks(list(np.arange(1.4, 3.8, .3)))
ax.set_yticks(list(np.arange(1.4, 3.8, .3)))
ax.set_xlim(1.3, 3.5)
ax.set_ylim(1.3, 3.5)
ax.tick_params(axis='x', length=4, width=1, labelsize=8)
ax.tick_params(axis='y', length=4, width=1, labelsize=8)
ax.set_xlabel('Original Promiscuity Values', fontsize=12)
ax.set_ylabel('Recalculated Promiscuity Values', fontsize=12)
plt.savefig(cwd / 'figure_2' / 'reproduction.pdf')
