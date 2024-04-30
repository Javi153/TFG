import aprox_hpnx as apx
import genetic_hpnx as ghpnx
import random
import math
import time

file = open("result3.txt", "w")
#filein = open("ejemplos.txt", "r")
res_apx = [0] * 200
res_gen = [0] * 200
time_apx = [0] * 200
time_gen = [0] * 200
for i in range(200):
    #str_seq = filein.readline()[:-1]
    #print('CASO %i: %s' % (i, str_seq))
    str_seq = ''.join([random.choice(['H', 'P', 'N', 'X']) for _ in range(random.randint(200, 500))])
    inicio = time.time()
    fin = 0
    _, res_apx[i] = apx.prot_fold(str_seq)
    fin = time.time()
    time_apx[i] = fin - inicio
    inicio = time.time()
    res_gen[i] = ghpnx.protein_fold(str_seq)
    fin = time.time()
    time_gen[i] = fin - inicio
#filein.close()
file.write("Resultado de aproximado: ")
file.write(' '.join(str(i) for i in res_apx))
file.write("\nTiempo de aproximado: ")
file.write(' '.join(str(i) for i in time_apx))
file.write("\nResultado de genético: ")
file.write(' '.join(str(i) for i in res_gen))
file.write("\nTiempo de genético: ")
file.write(' '.join(str(i) for i in time_gen))
file.close()