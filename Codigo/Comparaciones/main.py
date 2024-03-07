import aprox as apx
import genetic_hp as ghp
import random
import math

for i in range(100):
    r = random.randint(4, 50)
    str_seq = ''
    for _ in range(r):
        s = random.choice(['H', 'P'])
        str_seq += s
    print('CASO %i: %s' % (i, str_seq))
    apx.prot_fold(str_seq, 'C', math.sqrt)
    ghp.protein_fold(str_seq)