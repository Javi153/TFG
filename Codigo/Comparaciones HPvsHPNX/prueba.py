import random
import pandas as pd
import matplotlib.pyplot as plt

file = open("input2.txt", 'r')
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
        "20-50 aminoácidos": mean_scores[15:20] + mean_scores[0:5],
        "50-200 aminoácidos": mean_scores[20:25] + mean_scores[5:10],
        "200-500 aminoácidos": mean_scores[25:30] + mean_scores[10:15]
    },
    index=["Hib1_100", "Hib1_200", "Hib1_300", "Hib1_400", "Hib1_500", "Hib2_100", "Hib2_200", "Hib2_300", "Hib2_400", "Hib2_500"],
)
axis = dataframe.plot.bar(title = 'Algoritmos híbridos', ylabel = 'Energía', rot=0)
print(axis)
plt.show()"""


for i in range(5):
    print(sum([content_num[4*i + 2][j] <= content_num[4*i][j] for j in range(1000)]))