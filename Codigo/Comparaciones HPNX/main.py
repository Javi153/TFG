import aprox_hpnx as apx
import genetic_hpnx as ghpnx
import random
import math
import time

file = open("result6.txt", "w")
filein = open("input_large.txt", "r")
res_apx = [0] * 200
res_gen = [0] * 200
time_apx = [0] * 200
time_gen = [0] * 200
for i in range(200):
    str_seq = filein.readline()[:-1]
    print('CASO %i: %s' % (i, str_seq))
    #str_seq = ''.join([random.choice(['H', 'P', 'N', 'X']) for _ in range(random.randint(200, 500))])
    inicio = time.time()
    fin = 0
    _, res_apx[i] = apx.prot_fold(str_seq)
    fin = time.time()
    time_apx[i] = fin - inicio
    inicio = time.time()
    res_gen[i], time_gen[i] = ghpnx.protein_fold(str_seq)
    fin = time.time()
    time_gen[i] = [f - inicio for f in time_gen[i]]
filein.close()
file.write("Resultado de aproximado: ")
file.write(' '.join(str(i) for i in res_apx))
file.write("\nTiempo de aproximado: ")
file.write(' '.join(str(i) for i in time_apx))
for j in range(len(res_gen[0])):
    file.write("\nResultado de genético %i iteraciones: " % ((j+1) * 100))
    file.write(' '.join(str(i[j]) for i in res_gen))
    file.write("\nTiempo de genético %i iteraciones: "% ((j+1) * 100))
    file.write(' '.join(str(i[j]) for i in time_gen))
file.close()