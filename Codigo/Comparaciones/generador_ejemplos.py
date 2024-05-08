import random

file = open("input_very_large.txt", "w")
for i in range(200):
    prot = ''.join([random.choice(['H', 'P']) for _ in range(random.randint(500, 2000))]) + '\n'
    file.write(prot)

file.close()