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

mean_scores = []
mean_times = []
for i in range(len(content_num)):
    if i % 2 == 0:
        mean_scores.append(sum(content_num[i]) / len(content_num[i]))
    else:
        mean_times.append(sum(content_num[i]) / len(content_num[i]))

"""dataframe = pd.DataFrame(
    {
        "Aprox": mean_scores[:18:6],
        "Gen100": mean_scores[1:18:6],
        "Gen200": mean_scores[2:18:6],
        "Gen300": mean_scores[3:18:6],
        "Gen400": mean_scores[4:18:6],
        "Gen500": mean_scores[5:18:6]
    },
    index=["20-50 aminoácidos", "50-200 aminoácidos", "200-500 aminoácidos"],
)
axis = dataframe.plot.bar(title = 'Algoritmos para el modelo HP', ylabel = 'Energía', rot=0)
print(axis)
plt.show()

dataframe = pd.DataFrame(
    {
        "Aprox": mean_times[:18:6],
        "Gen100": mean_times[1:18:6],
        "Gen200": mean_times[2:18:6],
        "Gen300": mean_times[3:18:6],
        "Gen400": mean_times[4:18:6],
        "Gen500": mean_times[5:18:6]
    },
    index=["20-50 aminoácidos", "50-200 aminoácidos", "200-500 aminoácidos"],
)
axis = dataframe.plot.bar(title = 'Algoritmos para el modelo HP', ylabel = 'Tiempo (seg)', rot=0)
print(axis)
plt.show()

dataframe = pd.DataFrame(
    {
        "Aprox": mean_scores[18::6],
        "Hib100": mean_scores[19::6],
        "Hib200": mean_scores[20::6],
        "Hib300": mean_scores[21::6],
        "Hib400": mean_scores[22::6],
        "Hib500": mean_scores[23::6]
    },
    index=["20-50 aminoácidos", "50-200 aminoácidos", "200-500 aminoácidos"],
)
axis = dataframe.plot.bar(title = 'Algoritmo híbrido para el modelo HP', ylabel = 'Energía', rot=0)
print(axis)
plt.show()

dataframe = pd.DataFrame(
    {
        "Aprox": mean_times[:18:6],
        "Hib100": mean_times[1:18:6],
        "Hib200": mean_times[2:18:6],
        "Hib300": mean_times[3:18:6],
        "Hib400": mean_times[4:18:6],
        "Hib500": mean_times[5:18:6]
    },
    index=["20-50 aminoácidos", "50-200 aminoácidos", "200-500 aminoácidos"],
)
axis = dataframe.plot.bar(title = 'Algoritmo híbrido para el modelo HP', ylabel = 'Tiempo (seg)', rot=0)
print(axis)
plt.show()"""

aux = mean_scores[18:]
aux2 = mean_scores[:18]
for i in range(0, len(aux)-1):
    print("Comparacion %i: %f" % (i+1, aux[i+1] / aux[i]))
print("FLAG")

for i in range(0, len(aux)):
    print("Comparacion %i: %f" % (i+1, aux[i] / aux2[i]))

"""print(sum([content_num[34][i] <= content_num[70][i] for i in range(200)]))"""