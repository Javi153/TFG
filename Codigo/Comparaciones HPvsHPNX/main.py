import aprox_hpnx as hpnx_apx
import aprox as hp_apx
import random
import math
import time

file = open("result5.txt", "w")
#filein = open("ejemplos.txt", "r")
res_apx = [0] * 1000
res_gen = [0] * 1000
time_apx = [0] * 1000
time_gen = [0] * 1000
for i in range(1000):
    #str_seq = filein.readline()[:-1]
    #print('CASO %i: %s' % (i, str_seq))
    str_seq = ''.join([random.choice(['H', 'P', 'N', 'X']) for _ in range(random.randint(3000, 7000))])
    str_seq_hp = ''.join(['H' if str_seq[i] == 'H' else 'P' for i in range(len(str_seq))])
    inicio = time.time()
    fin = 0
    (ind_to_coord, coord_to_ind, coord), _ = hp_apx.prot_fold(str_seq_hp, 'C', hp_apx.f)
    res_apx[i] = hpnx_apx.fitness(str_seq, ind_to_coord, coord_to_ind)
    fin = time.time()
    time_apx[i] = fin - inicio
    inicio = time.time()
    _, res_gen[i] = hpnx_apx.prot_fold(str_seq)
    fin = time.time()
    time_gen[i] = fin - inicio
#filein.close()
file.write("Resultado de aproximado HP: ")
file.write(' '.join(str(i) for i in res_apx))
file.write("\nTiempo de aproximado HP: ")
file.write(' '.join(str(i) for i in time_apx))
file.write("\nResultado de aproximado HPNX: ")
file.write(' '.join(str(i) for i in res_gen))
file.write("\nTiempo de aproximado HPNX: ")
file.write(' '.join(str(i) for i in time_gen))
file.close()