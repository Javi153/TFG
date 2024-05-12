import random

file = open("input_small.txt", "w")
for i in range(1000):
    prot = ''.join([random.choice(['H', 'P', 'N', 'X']) for _ in range(random.randint(20, 50))]) + '\n'
    file.write(prot)
file.close()