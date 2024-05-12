import scikit_posthocs as sp
import numpy as np
from scipy import stats as st
import pandas as pd
import matplotlib.pyplot as plt

with open('input.txt', 'r') as infile:
    samples = [0] * 30
    for i in range(60):
        aux = infile.readline()[:-1]
        if i % 2 == 0:
            aux = aux.split(' ')
            samples[i // 2] = [int(x) for x in aux]
    samples = np.array(samples)
    for i in range(3):
        print('Caso %i' % i)
        pvalues = sp.posthoc_conover_friedman(np.concatenate((samples[5*i:5*i+5], samples[5*i+15:5*i+20])).T)
        for j in range(len(pvalues.values)):
            for num in range(len(pvalues[0])):
                pvalues[j][num] = "{:.3E}".format(pvalues[j][num])
        fig, ax = plt.subplots(1,1)
        ax.axis('off')
        pvalues.columns = ['Hib1_100', 'Hib1_200', 'Hib1_300', 'Hib1_400', 'Hib1_500', 'Hib2_100', 'Hib2_200', 'Hib2_300', 'Hib2_400', 'Hib2_500']
        pvalues.index = ['Hib1_100', 'Hib1_200', 'Hib1_300', 'Hib1_400', 'Hib1_500', 'Hib2_100', 'Hib2_200', 'Hib2_300', 'Hib2_400', 'Hib2_500']
        table = ax.table(cellText = pvalues.values, rowLabels = pvalues.index, colLabels = pvalues.columns, loc='center', cellLoc='center')
        table.auto_set_font_size(False)
        table.set_fontsize(12)
        plt.show()
        print(pvalues)
        stats, pvalue = st.friedmanchisquare(samples[5*i], samples[5*i+1], samples[5*i+2], samples[5*i+3], samples[5*i+4], samples[5*i+15], samples[5*i+16], samples[5*i+17], samples[5*i+18], samples[5*i+19])
        print(stats, pvalue)
    """for i in range(5):
        print('Caso %i' % i)
        pvalues = sp.posthoc_conover_friedman(samples[2*i:2*i+2].T)
        print(pvalues)"""
