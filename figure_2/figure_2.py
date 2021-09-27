##script for Figure 2A,B,C
import pandas as pd
import pylab as plt
import numpy as np
from kl import kl
from scipy.stats import spearmanr
from matplotlib import lines
import pathlib

path = pathlib.Path.cwd()
if path.stem == 'hla_promiscuity':
    cwd = path
else:
    cwd = list(path.parents)[::-1][path.parts.index('hla_promiscuity')]
    import sys
    sys.path.append(str(cwd))

iupac = 'ACDEFGHIKLMNPQRSTVWY'
cols = ['Reference IRI', 'Description', 'Date', 'Qualitative Measure', 'Allele Name', 'Method/Technique']
mhc_all = pd.read_csv(cwd / 'data' / 'mhc_ligand_full.csv', sep=',', skiprows=1, usecols=cols) ##http://www.iedb.org/database_export_v3.php
tcell_all = pd.read_csv(cwd / 'data' / 'tcell_full_v3.csv', sep=',', skiprows=1, usecols=cols) ##http://www.iedb.org/database_export_v3.php
tcell_all['Method/Technique'] = 'T cell'
mhc_all = mhc_all.append(tcell_all)
mhc_all['Reference IRI'] = mhc_all['Reference IRI'].apply(lambda x: x.split('/')[-1])
mhc_all = mhc_all.loc[mhc_all['Qualitative Measure'].str.contains('Positive')]
mhc_all = mhc_all.loc[mhc_all['Allele Name'].fillna('nan').str.contains('HLA-[A|B|C]...:..$')]
mhc_all['Description'] = mhc_all['Description'].str.upper()
mhc_all = mhc_all.loc[mhc_all['Description'].str.contains('^['+iupac+']+$')]
mhc_all = mhc_all.loc[(mhc_all['Description'].str.len() < 13) & (mhc_all['Description'].str.len() > 7)]
mhc_all['label'] = mhc_all['Allele Name'].apply(lambda x: x[4:].replace('*', '').replace(':', ''))


##peptides used by Manczinger et al.
mhc_df = pd.read_csv(cwd / 'data' / 'mhc_ligand.txt', sep='\t') ##https://www.nature.com/articles/s43018-021-00226-4
mhc_df['Description'] = mhc_df['Description'].str.upper()


##create dataframe for IEDB data used by Manczinger et al.
mhc_original = mhc_all.loc[mhc_all['Date'] <= 2018]
mhc_original = mhc_original.loc[mhc_original['Description'].isin(mhc_df['Description'].unique()) & mhc_original['Allele Name'].isin(mhc_df['Allele Name'].unique())]


original_reference_peptides_count = mhc_original.groupby(['label', 'Method/Technique'])['Description'].apply(lambda x: len(set(x))).to_dict()

fig1_data = pd.read_csv(cwd / 'data' / 'fig1_source_data.csv')  ##https://www.nature.com/articles/s43018-021-00226-4
original = {i: j for i, j in zip(fig1_data['Allele name'].values, fig1_data['Allele promiscuity'].values)}
original_sorted, sorted_names = zip(*sorted(list(zip(original.values(), original.keys()))))

original_mass_spec_counts = {}
original_T_cell_counts = {}
original_assay_counts = {}
for i in original_reference_peptides_count:
    if 'mass' in i[1]:
        original_mass_spec_counts[i[0]] = original_mass_spec_counts.get(i[0], 0) + original_reference_peptides_count[i]
    elif i[1] == 'T cell':
        original_T_cell_counts[i[0]] = original_T_cell_counts.get(i[0], 0) + original_reference_peptides_count[i]
    else:
        original_assay_counts[i[0]] = original_assay_counts.get(i[0], 0) + original_reference_peptides_count[i]

original_mass_spec_count_data = np.array([original_mass_spec_counts.get(i, 0) for i in sorted_names])
original_T_cell_count_data = np.array([original_T_cell_counts.get(i, 0) for i in sorted_names])
original_assay_count_data = np.array([original_assay_counts.get(i, 0) for i in sorted_names])

original_mass_spec_fraction = original_mass_spec_count_data / (original_mass_spec_count_data + original_T_cell_count_data + original_assay_count_data)
original_T_cell_count_fraction = original_T_cell_count_data / (original_mass_spec_count_data + original_T_cell_count_data + original_assay_count_data)
original_assay_count_fraction = original_assay_count_data / (original_mass_spec_count_data + original_T_cell_count_data + original_assay_count_data)


##calculate new pr values and the peptide fractions


##add pyke dataset, https://www.sciencedirect.com/science/article/pii/S1535947621000839
cols = ['peptide', 'allele', 'modification']
pyke = pd.read_csv(cwd / 'data' / '1-s2.0-S1535947621000839-mmc2-1.csv', sep=',', usecols=cols)
pyke = pyke.loc[pyke['modification'] == 'na']
pyke.drop('modification', axis=1, inplace=True)
pyke.rename(mapper={'allele': 'Allele Name', 'peptide': 'Description'}, axis=1, inplace=True)
pyke['Reference IRI'] = 'pyke'
pyke['Allele Name'] = pyke['Allele Name'].apply(lambda x: x[:5] + '*' + x[5:])
pyke['Method/Technique'] = 'mass spectrometry'
pyke['label'] = pyke['Allele Name'].apply(lambda x: x[4:].replace('*', '').replace(':', ''))
pyke['Description'] = pyke['Description'].str.upper()
pyke = pyke.loc[pyke['Description'].str.contains('^['+iupac+']+$')]
pyke = pyke.loc[(pyke['Description'].str.len() < 13) & (pyke['Description'].str.len() > 7)]

mhc_new = mhc_all.append(pyke)

new = {}
for allele in mhc_df['Allele Name'].unique():
    allele_df = mhc_new.loc[(mhc_new['Allele Name'] == allele)].drop_duplicates(subset=['Description'])
    if len(allele_df) > 400:
        new[allele] = kl(allele_df)


new_reference_peptides_count = mhc_new.groupby(['label', 'Method/Technique'])['Description'].apply(lambda x: len(set(x))).to_dict()

new_mass_spec_counts = {}
new_T_cell_counts = {}
new_assay_counts = {}
for i in new_reference_peptides_count:
    if 'mass' in i[1]:
        new_mass_spec_counts[i[0]] = new_mass_spec_counts.get(i[0], 0) + new_reference_peptides_count[i]
    elif i[1] == 'T cell':
        new_T_cell_counts[i[0]] = new_T_cell_counts.get(i[0], 0) + new_reference_peptides_count[i]
    else:
        new_assay_counts[i[0]] = new_assay_counts.get(i[0], 0) + new_reference_peptides_count[i]

new_mass_spec_count_data = np.array([new_mass_spec_counts.get(i, 0) for i in sorted_names])
new_T_cell_count_data = np.array([new_T_cell_counts.get(i, 0) for i in sorted_names])
new_assay_count_data = np.array([new_assay_counts.get(i, 0) for i in sorted_names])

new_mass_spec_fraction = new_mass_spec_count_data / (new_mass_spec_count_data + new_T_cell_count_data + new_assay_count_data)
new_T_cell_count_fraction = new_T_cell_count_data / (new_mass_spec_count_data + new_T_cell_count_data + new_assay_count_data)
new_assay_count_fraction = new_assay_count_data / (new_mass_spec_count_data + new_T_cell_count_data + new_assay_count_data)

fig = plt.figure()
gs = fig.add_gridspec(2, 1, height_ratios=[10, 10])
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[1, 0])

fig.subplots_adjust(top=1.00,
bottom=0.105,
left=0.075,
right=1.0,
hspace=0.09,
wspace=0.2)

original_values = np.array(original_sorted)
new_values = np.array([new['HLA-' + i[0] + '*' + i[1:-2] + ':' + i[-2:]][0] for i in sorted_names])
ax1.bar(np.arange(0, len(original_sorted)), original_values * original_mass_spec_fraction, width=.9, label='Mass Spec')
ax1.bar(np.arange(0, len(original_sorted)), original_values * original_assay_count_fraction, width=.9, bottom=original_values * original_mass_spec_fraction, label='HLA Binding')
ax1.bar(np.arange(0, len(original_sorted)), original_values * original_T_cell_count_fraction, width=.9, bottom=original_values * original_mass_spec_fraction + original_values * original_assay_count_fraction, label='T cell')
ax2.bar(np.arange(0, len(original_sorted)), new_values * new_mass_spec_fraction, width=.9, label='Mass Spec')
ax2.bar(np.arange(0, len(original_sorted)), new_values * new_assay_count_fraction, width=.9, bottom=new_values * new_mass_spec_fraction, label='HLA Binding')
ax2.bar(np.arange(0, len(original_sorted)), new_values * new_T_cell_count_fraction, width=.9, bottom=new_values * new_mass_spec_fraction + new_values * new_assay_count_fraction, label='T cell')


for ax in [ax1, ax2]:
    ax.spines['right'].set_visible(False)
    ax.spines['top'].set_visible(False)
    ax.spines['bottom'].set_visible(False)
    ax.spines['left'].set_linewidth(.8)
    ax.spines['left'].set_position(['outward', 2])
    ax.set_xticks([])
    ax.set_xlim(-.5, len(original_values) + 1)
    ax.tick_params(axis='y', length=3, width=.8, labelsize=6.5)
    ax.set_yticks(list(np.arange(0, 4, .5)))
    ax.spines['left'].set_bounds(0, 3.5)
    ax.set_ylim(0, max(original_values) + .2)
    ax.set_ylabel('Allele promiscuity', fontsize=10)
ax1.legend(frameon=False, ncol=3, loc=(.1, .88))


ax2.set_xticks(np.arange(0, len(sorted_names)))
ax2.set_xticklabels(sorted_names, rotation=90)
ax2.tick_params(axis='x', length=3, width=.8, labelsize=6.5, pad=0)
colors = ['#e5513b' if 'A' in i else '#3f5689' if 'B' in i else '#92d2c1' for i in sorted_names]

[t.set_color(colors[index]) for index, t in enumerate(ax2.xaxis.get_ticklabels())]

rect1 = lines.Line2D([], [], marker="s", markersize=4, linewidth=0, color='#e5513b')
rect2 = lines.Line2D([], [], marker="s", markersize=4, linewidth=0, color='#3f5689')
rect3 = lines.Line2D([], [], marker="s", markersize=4, linewidth=0, color='#92d2c1')
ax2.legend((rect1, rect2, rect3),
           ('HLA-A', 'HLA-B', 'HLA-C'),
           loc=(.4, -.25),
           frameon=False,
           borderpad=0,
           handletextpad=0,
           labelspacing=.3,
           fontsize=6,
           ncol=3)

plt.savefig(cwd / 'figure_2' / 'figure_2ab.pdf')

fig = plt.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(top=0.985,
bottom=0.1,
left=0.09,
right=0.975,
hspace=0.2,
wspace=0.2)

scatter = ax.scatter([original[i] for i in sorted_names], [new['HLA-' + i[0] + '*' + i[1:-2] + ':' + i[-2:]][0] for i in sorted_names],
                     edgecolor='k',)
ax.plot(np.arange(0, len(sorted_names) + 1), np.arange(0, len(sorted_names) + 1), color='k', alpha=.2, linestyle='--', linewidth=2)
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
ax.tick_params(length=4, width=1, labelsize=8)
ax.set_xlabel('Original Promiscuity Values', fontsize=12)
ax.set_ylabel('New Promiscuity Values', fontsize=12)
plt.savefig(cwd / 'figure_2' / 'figure_2c.pdf')

print(spearmanr([original[i] for i in original], [new['HLA-' + i[0] + '*' + i[1:-2] + ':' + i[-2:]][0] for i in original]))


##stem plot
mass_spec_delta = new_mass_spec_fraction - original_mass_spec_fraction
pr_delta = np.array([new['HLA-' + i[0] + '*' + i[1:-2] + ':' + i[-2:]][0] for i in sorted_names]) - np.array([original[i] for i in sorted_names])

fig = plt.figure()
ax = fig.add_subplot(111)
fig.subplots_adjust(top=1.0,
bottom=0.085,
left=0.09,
right=1.0,
hspace=0.2,
wspace=0.2)
ax.spines['right'].set_visible(False)
ax.spines['top'].set_visible(False)
ax.spines['left'].set_visible(False)
ax.spines['bottom'].set_visible(False)
ax.tick_params(width=0, labelsize=8, pad=-5)
ax.set_xlabel('New Minus Original Mass Spec Proportions', fontsize=12, labelpad=5)
ax.set_ylabel('New Minus Original Promiscuity Values', fontsize=12)
ax.set_xlim(min(mass_spec_delta) - .02, max(mass_spec_delta) + .02)
(markers, stemlines, baseline) = ax.stem(mass_spec_delta, pr_delta)
plt.setp(markers, markersize=4)
plt.setp(baseline, linestyle="-", color="grey", linewidth=1)
plt.setp(stemlines, linestyle="--", color="grey", linewidth=1)

plt.savefig(cwd / 'figure_2' / 'figure_2d.pdf')
