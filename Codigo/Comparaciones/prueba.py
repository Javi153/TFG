import random

file = open("input.txt", 'r')
content = file.readlines()
content_num = []
for line in content:
    line = line[:-1]
    line = line.split(' ')
    line = list(map(float, line))
    content_num.append(line)
    print(sum(line) / len(line))
print(sum(content_num[1]) + sum(content_num[3]))
aux = [content_num[0][i] / content_num[2][i] for i in range(len(content_num[0])) if content_num[2][i] != 0.0]
print(sum(aux) / len(aux))
file.close()