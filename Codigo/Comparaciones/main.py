import aprox as apx
import genetic_hp as ghp
import random
import math
import time

file = open("result.txt", "w")
res_apx = [0] * 100
res_gen = [0] * 100
time_apx = [0] * 100
time_gen = [0] * 100
for i in range(100):
    r = random.randint(500, 1000)
    str_seq = ''
    for _ in range(r):
        s = random.choice(['H', 'P'])
        str_seq += s
    print('CASO %i: %s' % (i, str_seq))
    inicio = time.time()
    fin = 0
    res_apx[i] = apx.prot_fold(str_seq, 'C', math.sqrt)
    fin = time.time()
    time_apx[i] = fin - inicio
    inicio = time.time()
    res_gen[i] = ghp.protein_fold(str_seq)
    fin = time.time()
    time_gen[i] = fin - inicio
file.write("Resultado de aproximado: ")
file.write(' '.join(str(i) for i in res_apx))
file.write("\nTiempo de aproximado: ")
file.write(' '.join(str(i) for i in time_apx))
file.write("\nResultado de genético: ")
file.write(' '.join(str(i) for i in res_gen))
file.write("\nTiempo de genético: ")
file.write(' '.join(str(i) for i in time_gen))
file.close()