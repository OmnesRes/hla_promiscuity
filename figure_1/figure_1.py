import pandas as pd
from scipy.stats import spearmanr
from kl import kl
import pylab as plt
import seaborn as sns

import pathlib
path = pathlib.Path.cwd()
if path.stem == 'hla_promiscuity':
    cwd = path
else:
    cwd = list(path.parents)[::-1][path.parts.index('hla_promiscuity')]
    import sys
    sys.path.append(str(cwd))

iupac = 'ACDEFGHIKLMNPQRSTVWY'
cols = ['Reference IRI', 'Description', 'Method/Technique', 'Qualitative Measure', 'Cell Type', 'Allele Name', 'Allele Evidence Code']
mhc_all = pd.read_csv(cwd / 'data' / 'mhc_ligand_full.csv', sep=',', skiprows=1, usecols=cols) ##http://www.iedb.org/database_export_v3.php
mhc_all['Reference IRI'] = mhc_all['Reference IRI'].apply(lambda x: x.split('/')[-1])
mhc_all = mhc_all.loc[mhc_all['Qualitative Measure'].str.contains('Positive')]
mhc_all = mhc_all.loc[mhc_all['Allele Name'].str.contains('HLA-[A|B|C]...:..$')]

##add pyke dataset, https://www.sciencedirect.com/science/article/pii/S1535947621000839
cols = ['peptide', 'allele', 'modification']
pyke = pd.read_csv(cwd / 'data' / '1-s2.0-S1535947621000839-mmc2-1.csv', sep=',', usecols=cols)
pyke = pyke.loc[pyke['modification'] == 'na']
pyke.drop('modification', axis=1, inplace=True)
pyke.rename(mapper={'allele': 'Allele Name', 'peptide': 'Description'}, axis=1, inplace=True)
pyke['Qualitative Measure'] = 'Positive'
pyke['Reference IRI'] = 'pyke'
pyke['Allele Name'] = pyke['Allele Name'].apply(lambda x: x[:5] + '*' + x[5:])
pyke['Method/Technique'] = 'mass spectrometry'

mhc_all = mhc_all.append(pyke)
mhc_all['Description'] = mhc_all['Description'].str.upper()
mhc_all = mhc_all.loc[mhc_all['Description'].str.contains('^['+iupac+']+$')]

idx = ~mhc_all['Method/Technique'].str.contains('mass')
mhc_all.loc[idx, 'Reference IRI'] = 'assay'
mhc_all = mhc_all.loc[(mhc_all['Description'].str.len() < 13) & (mhc_all['Description'].str.len() > 7)]


references = ['1036600', '1034798', '1032403', 'pyke', '1036766', 'assay', 'mass spec']

reference_dicts = []
for index, reference in enumerate(references):
    if index == 6: ##get kl for all mass spec data
        reference_df = mhc_all.loc[mhc_all['Reference IRI'] != 'assay']
    else:
        reference_df = mhc_all.loc[mhc_all['Reference IRI'] == reference]
    reference_dict = {}
    for allele in reference_df['Allele Name'].unique():
        allele_df = reference_df.loc[(reference_df['Allele Name'] == allele)].drop_duplicates(subset=['Description'])
        if len(allele_df) > 400:
            reference_dict[allele] = kl(allele_df)
    reference_dicts.append(reference_dict)

comparisons = {}
for index1, reference1 in enumerate(reference_dicts):
    for index2, reference2 in enumerate(reference_dicts):
        if index1 != index2:
            if len([i for i in reference1 if i in reference2]) >= 8:
                if tuple(sorted((index1, index2))) not in comparisons:
                    comparisons[tuple(sorted((index1, index2)))] = ''
                    print(index1, index2)
                    print(spearmanr([reference1[i][0] for i in reference1 if i in reference2], [reference2[i][0] for i in reference1 if i in reference2]),
                    len([i for i in reference1 if i in reference2]))


##scatter plots
labels = ['HLA Binding/Mass Spec', 'Sarkizova/Marcu', 'Pyke/Marcu', 'Sarkizova/Faridi', 'Sarkizova/Di Marco', 'Sarkizova/Pyke', 'Faridi/Pyke']
text = [['.02', '.92', '34'], ['.63', '.0002', '30'], ['.61', '.02', '14'], ['.29', '.28', '16'], ['.63', '.02', '14'], ['.61', '.004', '20'], ['-.07', '.87', '8']]
fig = plt.figure()
gs = fig.add_gridspec(3, 3)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])
ax4 = fig.add_subplot(gs[1, 0])
ax5 = fig.add_subplot(gs[1, 1])
ax6 = fig.add_subplot(gs[1, 2])
ax7 = fig.add_subplot(gs[2, 0])
ax8 = fig.add_subplot(gs[2, 1:])

fig.subplots_adjust(top=0.99,
bottom=0.08,
left=0.035,
right=1.0,
hspace=0.2,
wspace=0.34)


for index, (axis, data, label, text) in enumerate(zip([ax1, ax2, ax3, ax4, ax5, ax6, ax7], [[5, 6], [0, 4], [3, 4], [0, 1], [0, 2], [0, 3], [1, 3]], labels, text)):

    sns.regplot([reference_dicts[data[0]][i][0] for i in reference_dicts[data[0]] if i in reference_dicts[data[1]]],
                [reference_dicts[data[1]][i][0] for i in reference_dicts[data[0]] if i in reference_dicts[data[1]]],
                ax=axis,
                truncate=False,
                scatter_kws={'s': 5,
                             'linewidths': 0},
                line_kws={'lw': 1})
    axis.text(x=0.005, y=.9, s=r'$\rho$ =' + text[0] + r', p=' + text[1] + r', N=' + text[2], fontsize=5.5, transform=axis.transAxes)
    x_min = min([reference_dicts[data[0]][i][0] for i in reference_dicts[data[0]] if i in reference_dicts[data[1]]])
    x_max = max([reference_dicts[data[0]][i][0] for i in reference_dicts[data[0]] if i in reference_dicts[data[1]]])
    y_min = min([reference_dicts[data[1]][i][0] for i in reference_dicts[data[0]] if i in reference_dicts[data[1]]])
    y_max = max([reference_dicts[data[1]][i][0] for i in reference_dicts[data[0]] if i in reference_dicts[data[1]]])
    axis.set_xlim(x_min - (x_max - x_min) / 20, x_max + (x_max - x_min) / 20)
    axis.set_ylim(y_min - (y_max - y_min) / 20, y_max + (y_max - y_min) / 20)
    axis.set_xticks([x_min, x_max])
    axis.set_yticks([y_min, y_max])
    axis.set_xticklabels([round(i, 1) for i in [x_min, x_max]], fontsize=5.5)
    axis.set_yticklabels([round(i, 1) for i in [y_min, y_max]], fontsize=5.5)
    axis.set_xlabel(label.split('/')[0], fontsize=8, labelpad=-10)
    axis.set_ylabel(label.split('/')[1], fontsize=8, labelpad=-10)
    if "Marcu" in label:
        axis.yaxis.label.set_color("#bfbfbf")
    axis.tick_params(width=1, length=3)
    axis.spines['right'].set_visible(False)
    axis.spines['top'].set_visible(False)
    axis.spines['left'].set_bounds(y_min, y_max)
    axis.spines['bottom'].set_bounds(x_min, x_max)

##add monoallelic plot
references_to_use = ['1036600', '1034798', 'pyke']
reference_dicts_to_use = [reference_dicts[0]] + [reference_dicts[1]] + [reference_dicts[3]]


results_df = pd.DataFrame(data={'value': [ref[allele][0] for ref in reference_dicts_to_use for allele in ref],
                                'reference': [references_to_use[index] for index, ref in enumerate(reference_dicts_to_use) for allele in ref],
                                'allele': [allele for ref in reference_dicts_to_use for allele in ref]})

counts = results_df['allele'].value_counts()
alleles_to_use = [allele for allele, value in zip(counts.index, counts.values) if value == 3]


##calculate kl for all references

tcell_all = pd.read_csv(cwd / 'data' / 'tcell_full_v3.csv', sep=',', skiprows=1, usecols=['Reference IRI', 'Description', 'Method/Technique', 'Qualitative Measure', 'Cell Type', 'Allele Name', 'Allele Evidence Code']) ##http://www.iedb.org/database_export_v3.php
tcell_all['Method/Technique'] = 'T cell'
mhc_all = mhc_all.append(tcell_all)
mhc_all = mhc_all.loc[mhc_all['Qualitative Measure'].str.contains('Positive')]
mhc_all = mhc_all.loc[mhc_all['Allele Name'].fillna('nan').str.contains('HLA-[A|B|C]...:..$')]
mhc_all['Description'] = mhc_all['Description'].str.upper()
mhc_all = mhc_all.loc[mhc_all['Description'].str.contains('^['+iupac+']+$')]

all_kl = {}
for allele in alleles_to_use:
    allele_df = mhc_all.loc[mhc_all['Allele Name'] == allele].drop_duplicates(subset=['Description'])
    all_kl[allele] = kl(allele_df)

kl_df = pd.DataFrame(data={'value': [all_kl[i][0] for i in all_kl],
                                'reference': ['all' for i in all_kl],
                                'allele': [i for i in all_kl]})

kl_df['y'] = -.5

plotting_df = results_df.loc[results_df['allele'].isin(alleles_to_use)]
y_dict = {i: index for index, i in enumerate(references_to_use)}
plotting_df['y'] = plotting_df['reference'].apply(lambda x: y_dict[x])
plotting_df = kl_df.append(plotting_df)

##reference info
reference_labels = {'1036600': 'Sarkizova/\nB cell',
                    '1034798': 'Faridi/\nC1R cells',
                    'pyke': 'Pyke/\nK562',
                    }

import matplotlib.ticker as ticker
colors = ['#f94144ff', '#f3722cff', '#f8961eff', '#f9844aff', '#f9c74fff', '#90be6dff', '#43aa8bff', '#4d908eff', '#577590ff', '#277da1ff']


for index, allele in enumerate(alleles_to_use):
    ax8.plot(plotting_df.loc[plotting_df['allele'] == allele]['value'].values[0], plotting_df.loc[plotting_df['allele'] == allele]['y'].values[0], marker='|', zorder=1, ms=10, color=colors[index])
    ax8.plot(plotting_df.loc[plotting_df['allele'] == allele]['value'].values[1:], plotting_df.loc[plotting_df['allele'] == allele]['y'].values[1:], 'o', alpha=.40, zorder=1, ms=4, mew=0, color=colors[index])
    ax8.plot(plotting_df.loc[plotting_df['allele'] == allele]['value'].values, plotting_df.loc[plotting_df['allele'] == allele]['y'].values, linewidth=0.5, linestyle='--', zorder=-1, color=colors[index])

ax8.vlines([1.5, 2, 2.5], [-.25] * 3, [2.25] * 3, linestyles='dotted', alpha=.1, color='k', zorder=-1000)
ax8.spines['right'].set_visible(False)
ax8.spines['top'].set_visible(False)
ax8.spines['bottom'].set_visible(False)
ax8.spines['left'].set_visible(False)
ax8.set_yticks([-.75, 0, 1, 2])
ax8.set_yticklabels(['All data'] + [reference_labels[reference] for reference in references_to_use])
ax8.xaxis.set_label_position('top')
ax8.xaxis.set_major_locator(ticker.MultipleLocator(.5))
ax8.xaxis.set_minor_locator(ticker.FixedLocator([kl_df.loc[kl_df['allele'] == allele]['value'].values[0] for allele in plotting_df['allele'].unique()]))
ax8.xaxis.set_minor_formatter(ticker.FixedFormatter([i[4:].replace('*', '').replace(':', '') for i in plotting_df['allele'].unique()]))
ax8.tick_params(axis="x", which="minor", width=0, length=0, top=False, bottom=True, rotation=305, labelsize=8, pad=0)
ax8.tick_params(axis="x", which="major", width=0, length=0, top=True, labeltop=True, bottom=False, labelbottom=False, labelsize=8, pad=-8)
ax8.tick_params(axis="y", which="major", width=0, length=0, labelsize=8, pad=-5)
[t.set_color(colors[index]) for index, t in enumerate(ax8.xaxis.get_ticklabels(minor=True))]

plt.savefig(cwd / 'figure_1' / 'figure_1.pdf')