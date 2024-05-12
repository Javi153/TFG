import random
import pandas as pd
import matplotlib.pyplot as plt

file = open("input.txt", 'r')
content = file.readlines()
content_num = []
for line in content:
    line = line[:-1]
    line = line.split(' ')
    line = list(map(float, line))
    content_num.append(line)

"""mean_scores = []
mean_times = []
for i in range(len(content_num)):
    if i % 2 == 0:
        mean_scores.append(sum(content_num[i]) / len(content_num[i]))
    else:
        mean_times.append(sum(content_num[i]) / len(content_num[i]))


dataframe = pd.DataFrame(
    {
        "20-50 aminoácidos": mean_scores[0:6],
        "50-200 aminoácidos": mean_scores[6:12],
        "200-500 aminoácidos": mean_scores[12:18]
    },
    index=["Aprox", "Gen100", "Gen200", "Gen300", "Gen400", "Gen500"],
)
axis = dataframe.plot.bar(title = 'Algoritmos para el modelo HPNX', ylabel = 'Energía', rot=0)
print(axis)
plt.show()

dataframe = pd.DataFrame(
    {
        "20-50 aminoácidos": mean_scores[18:24],
        "50-200 aminoácidos": mean_scores[24:30],
        "200-500 aminoácidos": mean_scores[30:36]
    },
    index=["Aprox", "Hib100", "Hib200", "Hib300", "Hib400", "Hib500"],
)
axis = dataframe.plot.bar(title = 'Algoritmo híbrido para el modelo HPNX', ylabel = 'Energía', rot=0)
print(axis)
plt.show()"""

for i in range(5):
    print(sum([content_num[2*i + 26][j] <= content_num[2*i+62][j] for j in range(200)]))