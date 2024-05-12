import scikit_posthocs as sp
import numpy as np
from scipy import stats as st
import pandas as pd
import matplotlib.pyplot as plt

with open('input.txt', 'r') as infile:
    pd.set_option('display.precision', 3)
    samples = [0] * 36
    for i in range(72):
        aux = infile.readline()[:-1]
        if i % 2 == 0:
            aux = aux.split(' ')
            samples[i // 2] = [int(x) for x in aux]
    samples = np.array(samples)
    """for i in range(6):
        print('Caso %i' % i)
        pvalues = sp.posthoc_conover_friedman(samples[6*i:6*i+6].T)
        print(pvalues)
        pvalues = sp.posthoc_conover_friedman(samples[6*i+1:6*i+6].T)
        print(pvalues)
        stats, pvalue = st.friedmanchisquare(samples[6*i], samples[6*i+1], samples[6*i+2], samples[6*i+3], samples[6*i+4], samples[6*i+5])
        print(stats, pvalue)
        stats, pvalue = st.friedmanchisquare(samples[6*i+1], samples[6*i+2], samples[6*i+3], samples[6*i+4], samples[6*i+5])
        print(stats, pvalue)"""
    for i in range(3):
        pvalues = sp.posthoc_conover_friedman(np.concatenate((samples[6*i+1:6*i+6], samples[6*i+19:6*i+24])).T)
        for j in range(len(pvalues.values)):
            for num in range(len(pvalues[0])):
                pvalues[j][num] = "{:.3E}".format(pvalues[j][num])
        fig, ax = plt.subplots(1,1)
        ax.axis('off')
        pvalues.columns = ['Gen100', 'Gen200', 'Gen300', 'Gen400', 'Gen500', 'Hib100', 'Hib200', 'Hib300', 'Hib400', 'Hib500']
        pvalues.index = ['Gen100', 'Gen200', 'Gen300', 'Gen400', 'Gen500', 'Hib100', 'Hib200', 'Hib300', 'Hib400', 'Hib500']
        table = ax.table(cellText = pvalues.values, rowLabels = pvalues.index, colLabels = pvalues.columns, loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(12)
        plt.show()
        print(pvalues)
        stats, pvalue = st.friedmanchisquare(samples[6*i+1], samples[6*i+2], samples[6*i+3], samples[6*i+4], samples[6*i+5], samples[6*i+19], samples[6*i+20], samples[6*i+21], samples[6*i+22], samples[6*i+23])
        print(stats, pvalue)