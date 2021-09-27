import numpy as np
import pandas as pd
from scipy.stats import entropy

import pathlib
path = pathlib.Path.cwd()
if path.stem == 'hla_promiscuity':
    cwd = path
else:
    cwd = list(path.parents)[::-1][path.parts.index('hla_promiscuity')]
    import sys
    sys.path.append(str(cwd))


AA = ['A', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'K', 'L', 'M', 'N', 'P', 'Q',
       'R', 'S', 'T', 'V', 'W', 'Y']

stat_all = [pd.read_csv(cwd / 'data' / 'stat_all' / 'stat_0.txt', sep='\t'),
            pd.read_csv(cwd / 'data' / 'stat_all' / 'stat_1.txt', sep='\t'),
            pd.read_csv(cwd / 'data' / 'stat_all' / 'stat_2.txt', sep='\t'),
            pd.read_csv(cwd / 'data' / 'stat_all' / 'stat_3.txt', sep='\t'),
            pd.read_csv(cwd / 'data' / 'stat_all' / 'stat_4.txt', sep='\t')]

kl_all = [pd.read_csv(cwd / 'data' / 'kl' / 'kl_0.txt', sep='\t'),
            pd.read_csv(cwd / 'data' / 'kl' / 'kl_1.txt', sep='\t'),
            pd.read_csv(cwd / 'data' / 'kl' / 'kl_2.txt', sep='\t'),
            pd.read_csv(cwd / 'data' / 'kl' / 'kl_3.txt', sep='\t'),
            pd.read_csv(cwd / 'data' / 'kl' / 'kl_4.txt', sep='\t')]

##a function for calculating kl with the possibility of sampling
def kl(df, sampling=False, n=25, size=400):
    if sampling:
        n = n
    else:
        n = 1
    assay_kls = []
    for i in range(n):
        if sampling:
            indexes = np.random.permutation(np.arange(len(df)))[:size]
            subset_df = df.iloc[indexes, :]
        else:
            subset_df = df
        sizes = []
        all_entropies = []
        for index, size in enumerate(range(8, 13)):
            peptides = subset_df.loc[subset_df['Description'].str.len() == size]['Description'].apply(lambda x: list(x)).values
            if len(peptides) > 40:
                peptides = np.stack(peptides, axis=0)
                entropies = []
                for pos in range(size):
                    aa, counts = np.unique(peptides[:, pos], return_counts=True)
                    mask = [i in AA for i in aa]
                    aa = aa[mask]
                    counts = counts[mask]
                    if len(aa) < 20:
                        counts = counts + 0.0000001
                        freq = [counts[np.where(aa == i)[0][0]] if i in aa else 0.0000001 for i in AA]
                    else:
                        freq = counts
                    freq = np.array(freq) / sum(freq)
                    entr = entropy(freq, stat_all[index].iloc[:, pos].values)
                    entropies.append(entr)
                if [i for i in entropies if i >= .04] == []:
                    continue
                sizes.append(len(peptides))
                unnormalized = np.mean([i for i in entropies if i >= .04])
                nearest = int(round(len(peptides) / 20, 0)) * 20
                if nearest >= 50000:
                    reference = kl_all[index].iloc[50000 // 20 - 1, 1]
                else:
                    choice = nearest // 20 - 1
                    reference = kl_all[index].iloc[choice, 1]
                    if len(peptides) > nearest:
                        reference_2 = kl_all[index].iloc[choice + 1, 1]
                        reference = reference - np.abs(reference - reference_2) * (len(peptides) - nearest) / 20
                    elif len(peptides) < nearest:
                        reference_2 = kl_all[index].iloc[choice - 1, 1]
                        reference = reference + np.abs(reference - reference_2) * (nearest - len(peptides)) / 20
                    else:
                        reference = reference
                normalized = unnormalized - reference / 2
                all_entropies.append(normalized)
        if sizes != []:
            kl = sum([i * j for i, j in zip(sizes, all_entropies)]) / sum(sizes)
            assay_kls.append(1 / kl)
    if assay_kls != []:
        return [np.mean(assay_kls), np.std(assay_kls)]
    else:
        raise ValueError