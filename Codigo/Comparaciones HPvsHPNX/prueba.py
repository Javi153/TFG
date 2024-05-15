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
        "AproxHP": mean_scores[::2],
        "AproxHPNX": mean_scores[1::2],
    },
    index=["20-50\naminoácidos", "50-200\naminoácidos", "200-500\naminoácidos", "500-2000\naminoácidos", "2000-5000\naminoácidos"],
)
axis = dataframe.plot.bar(title = 'Algoritmos aproximados', ylabel = 'Energía', rot=0)
print(axis)
plt.show()"""

"""dataframe = pd.DataFrame(
    {
        "Hib1_100": mean_scores[:15:5],
        "Hib2_100": mean_scores[15::5],
        "Hib1_200": mean_scores[1:15:5],
        "Hib2_200": mean_scores[16::5],
        "Hib1_300": mean_scores[2:15:5],
        "Hib2_300": mean_scores[17::5],
        "Hib1_400": mean_scores[3:15:5],
        "Hib2_400": mean_scores[18::5],
        "Hib1_500": mean_scores[4:15:5],
        "Hib2_500": mean_scores[19::5]
    },
    index=["20-50 aminoácidos", "50-200 aminoácidos", "200-500 aminoácidos"],
)
axis = dataframe.plot.bar(title = 'Algoritmos híbridos', ylabel = 'Energía', rot=0)
print(axis)
plt.show()"""

for i in range(15):
    print(mean_scores[i+15] / mean_scores[i])

"""for i in range(5):
    print(sum([content_num[4*i + 2][j] <= content_num[4*i][j] for j in range(1000)]))"""