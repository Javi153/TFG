from enum import Enum
import numpy as np
import random
import aprox as apx
import math
import proteins as pts
from proteins import Directions
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d

#Directions = Enum('Directions', ['U', 'D', 'L', 'R', 'F', 'B'])

def next_dir(d: pts.Directions, c: tuple[int, int, int]):
    if d == pts.Directions.U:
        return (c[0], c[1] + 1, c[2])
    elif d == pts.Directions.D:
        return (c[0], c[1] - 1, c[2])
    elif d == pts.Directions.L:
        return (c[0] - 1, c[1], c[2])
    elif d == pts.Directions.R:
        return (c[0] + 1, c[1], c[2])
    elif d == pts.Directions.F:
        return (c[0], c[1], c[2] + 1)
    elif d == pts.Directions.B:
        return (c[0], c[1], c[2] - 1)

def neighbours(c: tuple[int, int, int]) -> list[pts.Directions]:
    res = []
    for d in pts.Directions:
        res.append(next_dir(d, c))
    return res

def free_neighbours(c: tuple[int, int, int], coor_dic: dict) -> list[pts.Directions]:
    res = []
    for d in pts.Directions:
        if next_dir(d, c) not in coor_dic:
            res.append(d)
    return res

def initial_population(p: str, cant: int) -> list[tuple[list[pts.Directions], dict, dict]]:
    res = []
    repetition = -1
    i_repeat = -1
    for _ in range(cant):
        dir = [0] * (len(p) - 1)
        dir[0] = random.choice(list(pts.Directions))
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

def fitness(p: str, ind_to_coor: dict, coor_to_ind: dict) -> int:
    energy = 0
    for i in range(len(p)):
        for n in neighbours(ind_to_coor[i]):
            if coor_to_ind.get(n) != None and coor_to_ind.get(n) > i and coor_to_ind.get(n) - i > 1 and p[i] == 'H' and p[coor_to_ind.get(n)] == 'H':
                energy += 1
    return energy

def mutation(p: list[pts.Directions], ind_to_coor: dict, coor_to_ind: dict) -> tuple[list[pts.Directions], dict, dict]:
    mut_point = random.randint(1, len(p) - 2)
    res = p[0:mut_point]
    res.append(random.choice(list(pts.Directions)))
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

def cross(p1: list[pts.Directions], p2: list[pts.Directions]) -> tuple[tuple[list[pts.Directions], dict, dict], tuple[list[pts.Directions], dict, dict]]:
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

def genetic(N: int, it: int, p: str):
    population = None
    while population == None:
        population = initial_population(p, N-1)
    res_aprox = apx.algorithmC(apx.format_seq(p), apx.f)
    ind_to_coor_aprox, coor_to_ind_aprox, _ = apx.prot_coord(res_aprox, in3D = True)
    population.append((res_aprox, ind_to_coor_aprox, coor_to_ind_aprox))
    best_ins = []
    best_sol = 0
    for _ in range(it):
        new_pop = []
        pop_fitness = [fitness(p, ind, coor) for (_, ind, coor) in population]
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
            if fitness(p, to_cross[i][1], to_cross[i][2]) < fitness(p, child1[1], child1[2]):
                new_pop.append(child1)
            else:
                new_pop.append(to_cross[i])
            if fitness(p, to_cross[i+1][1], to_cross[i+1][2]) < fitness(p, child2[1], child2[2]):
                new_pop.append(child2)
            else:
                new_pop.append(to_cross[i+1])
        population = new_pop
    pop_fitness = [fitness(p, ind, coor) for (_, ind, coor) in population]
    aux_sol = max(pop_fitness)
    if aux_sol > best_sol:
        best_sol = aux_sol
        best_ins = population[np.argmax(pop_fitness)]
    return best_ins, best_sol

def protein_fold(p: str):
    ins, fit = genetic(100, 150, p)
    print('El valor obtenido por el algoritmo genético es', fit)
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

#str_seq = 'HPPPHPPPHPPPHHPPPHPPPHPPPHPPPHPPPHPPPHPHPPPHPPPHPPPHHPPPHPPPHPPPHPPPHPPPHPPPHPHPPPHPPPHPPPHHPPPHPPPHPPPHPPPHPPPHPPPHPHPPPHPPPHPPPHHPPPHPPPHPPPHPPPHPPPHPPPHP'
#protein_fold(str_seq)

#str_seq = 'HPPHPHHPHHPPHPHPHHHPHPHPPHPHPHHPHPPHHHPHPHHPHPPPPHHHPPHHPHPHHPHHPPHPHPPHPPHHHPPHPPHPPHPHHPHPHPPHPHPHHPHPPHHHPHPHHPHPPPPHHHPPHHPHPHHPHH'
#protein_fold(str_seq)

#str_seq = 'HPPHHHPHPPPHHPHPPPHPHHHPPHHPHPHPHPHHHPPHPHPPHPHPHHHPHPHPHHHPHPPPHHPPHPPPHPPHPHPPPPHHPHPPHPHPHHPPHPHPPHPHHPPHHPHPHHPHPHPPHHHPHP'
#protein_fold(str_seq)