from enum import Enum
import numpy as np
import random
import math
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import subprocess
from subprocess import PIPE
import os

Directions = Enum('Directions', ['U', 'D', 'L', 'R', 'F', 'B'])

NO_GRUPOS = 5

def next_dir(d: Directions, c: tuple[int, int, int]):
    if d == Directions.U:
        return (c[0], c[1] + 1, c[2])
    elif d == Directions.D:
        return (c[0], c[1] - 1, c[2])
    elif d == Directions.L:
        return (c[0] - 1, c[1], c[2])
    elif d == Directions.R:
        return (c[0] + 1, c[1], c[2])
    elif d == Directions.F:
        return (c[0], c[1], c[2] + 1)
    elif d == Directions.B:
        return (c[0], c[1], c[2] - 1)

def neighbours(c: tuple[int, int, int]) -> list[Directions]:
    res = []
    for d in Directions:
        res.append(next_dir(d, c))
    return res

def free_neighbours(c: tuple[int, int, int], coor_dic: dict) -> list[Directions]:
    res = []
    for d in Directions:
        if next_dir(d, c) not in coor_dic:
            res.append(d)
    return res

def translation(p: str) -> str:
    result = ''
    for a in p:
        if a in ['ALA', 'VAL']:
            result += 'h'
        elif a in ['GLY', 'ILE', 'LEU', 'MET', 'PHE', 'PRO', 'TRP']:
            result += 'H'
        elif a in ['ARG', 'HIS', 'LYS']:
            result += 'P'
        elif a in ['ASP', 'GLU']:
            result += 'N'
        else:
            result += 'X'
    return result


def initial_population(p: str, cant: int) -> list[tuple[list[Directions], dict, dict]]:
    res = []
    repetition = -1
    i_repeat = -1
    for _ in range(cant):
        dir = [0] * (len(p) - 1)
        dir[0] = random.choice(list(Directions))
        next = next_dir(dir[0], (0,0,0))
        ind_to_coor = {0 : (0,0,0), 1 : next}
        coor_to_ind = {(0,0,0) : 0, next : 1}
        i = 1
        while i < len(p) - 1:
            neigh = free_neighbours(next, coor_to_ind)
            if not neigh:
                coor_to_ind.pop(next)
                ind_to_coor.pop(i)
                i -= 1
                next = ind_to_coor[i]
                if i_repeat != i:
                    i_repeat = i
                    repetition = 0
                else:
                    repetition += 1
                    if repetition > 10:
                        return None
            else:
                dir[i] = random.choice(neigh)
                next = next_dir(dir[i], next)
                i += 1
                coor_to_ind[next] = i
                ind_to_coor[i] = next
        res.append((dir, ind_to_coor, coor_to_ind))
    return res

def amino_to_ind(c: str) -> int:
    if c == 'h':
        return 0
    elif c == 'H':
        return 1
    elif c == 'P':
        return 2
    elif c == 'N':
        return 3
    else:
        return 4

def get_energy(p: str, i: int, j: int, table: list[list[int]]) -> int:
    first = amino_to_ind(p[i])
    second = amino_to_ind(p[j])
    if second > first:
        first, second = second, first
    return table[first][second]

def fitness(p: str, ind_to_coor: dict, coor_to_ind: dict, table: list[list[int]]) -> int:
    energy = 0
    for i in range(len(p)):
        for n in neighbours(ind_to_coor[i]):
            if coor_to_ind.get(n) != None and coor_to_ind.get(n) > i and coor_to_ind.get(n) - i > 1 and p[i] == 'H' and p[coor_to_ind.get(n)] == 'H':
                energy += get_energy(p, i, coor_to_ind.get(n), table)
    return energy

def mutation(p: list[Directions], ind_to_coor: dict, coor_to_ind: dict) -> tuple[list[Directions], dict, dict]:
    mut_point = random.randint(1, len(p) - 2)
    res = p[0:mut_point]
    res.append(random.choice(list(Directions)))
    res = res + p[mut_point + 1:]
    ind_to_coor = {0 : (0,0,0)}
    coor_to_ind = {(0,0,0) : 0}
    i = 0
    repetition = -1
    i_repeat = -1
    while i < len(p):
        aux_next = next_dir(res[i], ind_to_coor[i])
        if coor_to_ind.get(aux_next) != None:
            if not free_neighbours(ind_to_coor[i], coor_to_ind):
                coor_to_ind.pop(ind_to_coor[i])
                ind_to_coor.pop(i)
                i -= 1
            if i_repeat != i:
                i_repeat = i
                repetition = 0
            else:
                repetition += 1
                if repetition > 10:
                    return None
            res[i] = random.choice(free_neighbours(ind_to_coor[i], coor_to_ind))
        else:
            coor_to_ind[aux_next] = i + 1
            ind_to_coor[i + 1] = aux_next
            i += 1
    return (p, ind_to_coor, coor_to_ind)

def cross(p1: list[Directions], p2: list[Directions]) -> tuple[tuple[list[Directions], dict, dict], tuple[list[Directions], dict, dict]]:
    cross_point = random.randint(1, len(p1) - 2)
    child1 = [d for d in p1[0:cross_point]] + [d for d in p2[cross_point:]]
    child2 = [d for d in p2[0:cross_point]] + [d for d in p1[cross_point:]]
    ind_to_coor1 = {0 : (0,0,0)}
    ind_to_coor2 = {0 : (0,0,0)}
    coor_to_ind1 = {(0,0,0) : 0}
    coor_to_ind2 = {(0,0,0) : 0}
    i = 0
    repetition = -1
    i_repeat = -1
    while i < len(p1):
        aux_next = next_dir(child1[i], ind_to_coor1[i])
        if coor_to_ind1.get(aux_next) != None:
            if not free_neighbours(ind_to_coor1[i], coor_to_ind1):
                coor_to_ind1.pop(ind_to_coor1[i])
                ind_to_coor1.pop(i)
                i -= 1
            if i_repeat != i:
                i_repeat = i
                repetition = 0
            else:
                repetition += 1
                if repetition > 10:
                    return None, None
            child1[i] = random.choice(free_neighbours(ind_to_coor1[i], coor_to_ind1))
        else:
            coor_to_ind1[aux_next] = i + 1
            ind_to_coor1[i + 1] = aux_next
            i += 1
    i = 0
    while i < len(p1):
        aux_next = next_dir(child2[i], ind_to_coor2[i])
        if coor_to_ind2.get(aux_next) != None:
            if not free_neighbours(ind_to_coor2[i], coor_to_ind2):
                coor_to_ind2.pop(ind_to_coor2[i])
                ind_to_coor2.pop(i)
                i -= 1
            if i_repeat != i:
                i_repeat = i
                repetition = 0
            else:
                repetition += 1
                if repetition > 10:
                    return None, None
            child2[i] = random.choice(free_neighbours(ind_to_coor2[i], coor_to_ind2))
        else:
            coor_to_ind2[aux_next] = i + 1
            ind_to_coor2[i + 1] = aux_next
            i += 1
    return (child1, ind_to_coor1, coor_to_ind1), (child2, ind_to_coor2, coor_to_ind2)

def genetic(N: int, it: int, p: str, table: list[list[int]]):
    population = None
    while population == None:
        population = initial_population(p, N)
    best_ins = []
    best_sol = 0
    for _ in range(it):
        new_pop = []
        pop_fitness = [fitness(p, ind, coor, table) for (_, ind, coor) in population]
        aux_sol = max(pop_fitness)
        if aux_sol > best_sol:
            best_sol = aux_sol
            best_ins = population[np.argmax(pop_fitness)]
        to_cross = []
        for i in range(N // 2):
            r1 = random.randint(0, len(population) - 1)
            parent1 = population.pop(r1)
            r2 = random.randint(0, len(population) - 1)
            parent2 = population.pop(r2)
            if pop_fitness[r1] > pop_fitness[r2]:
                to_cross.append(parent1)
                new_pop.append(parent2)
            else:
                to_cross.append(parent2)
                new_pop.append(parent1)
        for i in range(0, len(to_cross), 2):
            if i + 1 >= len(to_cross):
                new_pop.append(to_cross[i])
                break
            r = random.random()
            if r <= 0.7:
                child1, child2 = None, None
                while child1 == None and child2 == None:
                    child1, child2 = cross(to_cross[i][0], to_cross[i + 1][0])
            else:
                child1, child2 = to_cross[i], to_cross[i + 1]
            r = random.random()
            if r <= 0.1:
                aux = None
                while aux == None:
                    aux = mutation(child1[0], child1[1], child1[2])
                child1 = aux
            r = random.random()
            if r <= 0.1:
                aux = None
                while aux == None:
                    aux = mutation(child2[0], child2[1], child2[2])
                child2 = aux
            if fitness(p, to_cross[i][1], to_cross[i][2], table) < fitness(p, child1[1], child1[2], table):
                new_pop.append(child1)
            else:
                new_pop.append(to_cross[i])
            if fitness(p, to_cross[i+1][1], to_cross[i+1][2], table) < fitness(p, child2[1], child2[2], table):
                new_pop.append(child2)
            else:
                new_pop.append(to_cross[i+1])
        population = new_pop
    pop_fitness = [fitness(p, ind, coor, table) for (_, ind, coor) in population]
    aux_sol = max(pop_fitness)
    if aux_sol > best_sol:
        best_sol = aux_sol
        best_ins = population[np.argmax(pop_fitness)]
    return best_ins, best_sol

def initial_tables(N):
    tables = [0] * N
    for i in range(N):
        tab = []
        for j in range(NO_GRUPOS):
            tab.append([random.randint(-10, 10) for _ in range(NO_GRUPOS - j)])
        tables[i] = tab
    return tables

def compare(ins, media: float, archivo) -> float:
    subprocess.run("cp %s ../LGA_package_src/MOL2/prueba.pdb" % archivo, shell = True, executable="/bin/bash", stdout=PIPE, stderr=PIPE)
    f = open("../LGA_package_src/MOL2/prueba.pdb", "a")
    i = 0
    f.write("\nMOLECULE genetico\n")
    for c in ins[2]:
        f.write("ATOM  ")
        f.write("%5d " % (i+1))
        f.write("CA   ")
        f.write("%s A%4d    " % ("ABC", i+1))
        f.write("%8.3f" % (c[0] * media))
        f.write("%8.3f" % (c[1] * media))
        f.write("%8.3f  1.00 83.00           C\n" % (c[2] * media))
        i += 1
    f.write("END\n")
    f.close()
    p1 = subprocess.run("""cd ../LGA_package_src/; ulimit -s unlimited; ./lga -3 -gdt prueba.pdb | grep "GDT PERCENT_AT" | awk '{ V=($4+$6+$10+$18)/4.0; printf \"%7.3f\",V; }'""", shell = True, stdout = PIPE)
    return float(p1.stdout)


def protein_fold(N):
    #leer proteinas, supongamos lp
    archivos = os.listdir(".")
    archivos.remove('supergenetic.py')
    tables = initial_tables(N)
    for archivo in archivos:
        coords = []
        p = []
        with open(archivo, "r") as file:
            while True:
                line = file.readline()
                if not line:
                    break
                line = line.split()
                if line[0] == 'ATOM' and line[2] == 'CA':
                    coords.append((float(line[6]), float(line[7]), float(line[8])))
                    p.append(line[3])
        media = sum([abs(coords[i][0] - coords[i-1][0]) + abs(coords[i][1] - coords[i-1][1]) + abs(coords[i][2] - coords[i-1][2]) for i in range(1, len(coords))]) / (len(coords) - 1)
        for i in range(N):
            print(tables[i])
            ins, fit = genetic(100, 500, translation(p), tables[i])
            print('El valor obtenido por el algoritmo genético con la tabla %d es %s' % (i, fit))
            print('El valor de similitud para la proteína %s es de %7.3f%%' % (archivo, compare(ins, media, archivo)))
    """color = []
    for ch in p:
        if ch == 'H':
            color.append('black')
        else:
            color.append('white')
    coord = [c for c in ins[2]]
    coord_x = [c[0] for c in coord]
    coord_y = [c[1] for c in coord] 
    coord_z = [c[2] for c in coord]
    ax = plt.axes(projection ='3d')
    ax.plot3D(coord_x, coord_y, coord_z, c = 'black', zorder = 1)
    ax.scatter3D(coord_x, coord_y, coord_z, c = color, zorder = 2)
    plt.axis('equal')
    plt.show()"""

protein_fold(1)

#str_seq = 'HPPPHPPPHPPPHHPPPHPPPHPPPHPPPHPPPHPPPHPHPPPHPPPHPPPHHPPPHPPPHPPPHPPPHPPPHPPPHPHPPPHPPPHPPPHHPPPHPPPHPPPHPPPHPPPHPPPHPHPPPHPPPHPPPHHPPPHPPPHPPPHPPPHPPPHPPPHP'
#protein_fold(str_seq)

#str_seq = 'HPPHPHHPHHPPHPHPHHHPHPHPPHPHPHHPHPPHHHPHPHHPHPPPPHHHPPHHPHPHHPHHPPHPHPPHPPHHHPPHPPHPPHPHHPHPHPPHPHPHHPHPPHHHPHPHHPHPPPPHHHPPHHPHPHHPHH'
#protein_fold(str_seq)

#str_seq = 'HPPHHHPHPPPHHPHPPPHPHHHPPHHPHPHPHPHHHPPHPHPPHPHPHHHPHPHPHHHPHPPPHHPPHPPPHPPHPHPPPPHHPHPPHPHPHHPPHPHPPHPHHPPHHPHPHHPHPHPPHHHPHP'
#protein_fold(str_seq)