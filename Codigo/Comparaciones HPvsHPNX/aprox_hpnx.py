import hpnx_proteins as hpnxpts
import proteins as pts
from matplotlib import pyplot as plt
from mpl_toolkits import mplot3d
import random

def Mxy(p1: hpnxpts.prot, p2: hpnxpts.prot, getmin = True) -> int:
    if getmin:
        return min(p1.Nx(), p2.Ny())
    else:
        return max(p1.Nx(), p2.Ny())
    
def Myx(p1: hpnxpts.prot, p2: hpnxpts.prot, getmin = True) -> int:
    return Mxy(p2, p1, getmin)
    
def cond(a, b, A, B) -> bool:
    return a < b or (a == b and A < B)

def subroutine1(p : hpnxpts.prot):
    if len(p.getBlocks()) < 3: 
        return p, hpnxpts.prot([], False), False
    reverse = False
    P1 = []
    P2 = []
    pb = hpnxpts.prot_block(hpnxpts.Block_type.SEP)
    pb.add_amino(0)
    B1 = hpnxpts.prot(p.getBlocks()[0:2] + [pb])
    sep = p.getBlocks()[2]
    B2 = hpnxpts.prot([pb] + p.getBlocks()[3:])
    m1 = Mxy(B1, B2)
    m2 = Myx(B1, B2)
    M1 = Mxy(B1, B2, False)
    M2 = Myx(B1, B2, False)
    if m1 > m2:
        P1 = B2
        P2 = B1
        e = m1
        E = M1
        reverse = True
    else:
        P1 = B1
        P2 = B2
        e = m2
        E = M2
        reverse = False
    for i in range(3, (len(p.getBlocks()) + 1) // 2):
        B1 = hpnxpts.prot(p.getBlocks()[0:2*(i-1)] + [pb])
        B2 = hpnxpts.prot([pb] + p.getBlocks()[(2*i-1):])
        m1 = Mxy(B1, B2)
        m2 = Myx(B1, B2)
        M1 = Mxy(B1, B2, False)
        M2 = Myx(B1, B2, False)
        if cond(e, m1, E, M1):
            P1 = B2
            P2 = B1
            e = m1
            E = M1
            reverse = True
            sep = p.getBlocks()[2*(i-1)]
        if cond(e, m2, E, M2):
            P1 = B1
            P2 = B2
            e = m2
            E = M2
            reverse = False
            sep = p.getBlocks()[2*(i-1)]
    return P1, sep, P2, reverse

def format_seq(str_seq: str) -> list[int]:
    next = ['P', 'N', 'X']
    cant = 0
    res = []
    for ch in str_seq:
        if ch in next:
            cant += 1
        else:
            res.append(cant)
            if next != ['H']:
                next = ['H']
            else:
                next = ['P', 'N', 'X']
            cant = 1
    if cant > 0:
        res.append(cant)
    if len(res) % 2 == 0:
        res.append(0)
    return res

def aproxHPNX(p: str):
    paux = format_seq(p)
    dir = []
    if sum(paux) <= 3 or len(paux) == 1:
        for _ in range(sum(paux) - 1):
            dir.append(pts.Directions.R)
        return dir
    paux = hpnxpts.prot(paux, False)
    pb = hpnxpts.prot_block(hpnxpts.Block_type.SEP)
    pb.add_amino(0)
    for _ in range(0, paux.getBlocks()[0].getSize()):
        dir.append(pts.Directions.D)
    P1, sep, P2, reverse = subroutine1(hpnxpts.prot([pb] + paux.getBlocks()[1:-1] + [pb]))
    if len(P1.getBlocks()) <= 1 or len(P2.getBlocks()) <= 1:
        return [pts.Directions.R for _ in range(len(p) - 1)]
    P1.setType(hpnxpts.Block_type.Y_BLOCK)
    P2.setType(hpnxpts.Block_type.X_BLOCK)
    if reverse:
        fold1 = P2.centrals()
        fold2 = P1.centrals()
        acum = paux.getBlocks()[0].getSize() + sum([bl.getSize() for bl in P2.getBlocks()]) + sep.getSize()
    else:
        fold1 = P1.centrals()
        fold2 = P2.centrals()
        acum = paux.getBlocks()[0].getSize() + sum([bl.getSize() for bl in P1.getBlocks()]) + sep.getSize()
    fold1 = list(map(lambda num: num + paux.getBlocks()[0].getSize(), fold1))
    fold2 = list(map(lambda num: num + acum, fold2))
    if fold1[0] != paux.getBlocks()[0].getSize():
        for _ in range(fold1[0] - paux.getBlocks()[0].getSize()):
            dir.append(pts.Directions.D)
    left_central = [pts.Directions.D for _ in range(2 * len(fold1) - 1)]
    right_central = [pts.Directions.U for _ in range(2 * len(fold2) - 1)]
    left_laterals_ind = []
    right_laterals_ind = []
    if len(fold1) <= len(fold2):
        for i in range(0, len(fold1) - 1):
            if (p[fold1[i] + 1] == 'P' and p[fold2[len(fold1) - i - 2] + 1] == 'P') or (p[fold1[i] + 1] == 'N' and p[fold2[len(fold1) - i - 2] + 1] == 'N'):
                left_central[2 * i] = pts.Directions.F
                right_central[len(left_central) - 1 - 2 * (i + 1)] = pts.Directions.B
    else:
        for i in range(0, len(fold2) - 1):
            if (p[fold1[len(fold1) - len(fold2) + i] + 1] == 'P' and p[fold2[len(fold2) - i - 2] + 1] == 'P') or (p[fold1[len(fold1) - len(fold2) + i] + 1] == 'N' and p[fold2[len(fold2) - i - 2] + 1] == 'N'):
                left_central[len(left_central) - len(right_central) + 2 * i] = pts.Directions.F
                right_central[-2 * i - 3] = pts.Directions.B
    for i in range(1, len(fold1)):
        left_laterals_ind.append(list(range(fold1[i - 1] + 2, fold1[i])))
    for i in range(1, len(fold2)):
        right_laterals_ind.append(list(range(fold2[i - 1] + 2, fold2[i])))
    """left_laterals_mod = []
    right_laterals_mod = []
    for l in left_laterals_ind:
        ind = 0
        plus_links = 0
        minus_links = 0
        while ind < len(l) - 1 - ind:
            if (p[l[ind]] == 'P' and p[l[len(l) - 1 - ind]] == 'N') or (p[l[ind]] == 'N' and p[l[len(l) - 1 - ind]] == 'P'):
                minus_links += 1
            elif (p[l[ind]] == 'P' and p[l[len(l) - 1 - ind]] == 'P') or (p[l[ind]] == 'N' and p[l[len(l) - 1 - ind]] == 'N'):
                plus_links += 1
            ind += 1
        if minus_links < plus_links:
            left_laterals_mod.append(True)
        else:
            left_laterals_mod.append(False)
    for l in right_laterals_ind:
        ind = 0
        plus_links = 0
        minus_links = 0
        while ind < len(l) - 1 - ind:
            if (p[l[ind]] == 'P' and p[l[len(l) - 1 - ind]] == 'N') or (p[l[ind]] == 'N' and p[l[len(l) - 1 - ind]] == 'P'):
                minus_links += 1
            elif (p[l[ind]] == 'P' and p[l[len(l) - 1 - ind]] == 'P') or (p[l[ind]] == 'N' and p[l[len(l) - 1 - ind]] == 'N'):
                plus_links += 1
            ind += 1
        if minus_links < plus_links:
            right_laterals_mod.append(True)
        else:
            right_laterals_mod.append(False)"""
    for i in range(0, len(fold1) - 1):
        dir.append(left_central[2 * i])
        l = left_laterals_ind[i]
        #if left_laterals_mod[i]:
        if len(l) > 3:
            dir.append(pts.Directions.F)
            for _ in range(len(l) // 2 - 1):
                dir.append(pts.Directions.L)
            dir.append(pts.Directions.B)
            dir.append(left_central[2 * i + 1])
            for _ in range(len(l) // 2 - 1):
                dir.append(pts.Directions.R)
        else:
            for _ in range(len(l) // 2):
                dir.append(pts.Directions.L)
            dir.append(left_central[2 * i + 1])
            for _ in range(len(l) // 2):
                dir.append(pts.Directions.R)
    for _ in range(sep.getSize() // 2):
        dir.append(pts.Directions.D)
    dir.append(pts.Directions.R)
    for _ in range(sep.getSize() // 2):
        dir.append(pts.Directions.U)
    for i in range(0, len(fold2) - 1):
        dir.append(right_central[2 * i])
        l = right_laterals_ind[i]
        #if right_laterals_mod[i]:
        if len(l) > 3:
            dir.append(pts.Directions.B)
            for _ in range(len(l) // 2 - 1):
                dir.append(pts.Directions.R)
            dir.append(pts.Directions.F)
            dir.append(right_central[2 * i + 1])
            for _ in range(len(l) // 2 - 1):
                dir.append(pts.Directions.L)
        else:
            for _ in range(len(l) // 2):
                dir.append(pts.Directions.R)
            dir.append(right_central[2 * i + 1])
            for _ in range(len(l) // 2):
                dir.append(pts.Directions.L)
    if fold2[-1] + paux.getBlocks()[-1].getSize() != len(p) - 1:
        for _ in range(len(p) - 1 - fold2[-1] - paux.getBlocks()[-1].getSize()):
            dir.append(pts.Directions.U)
    for _ in range(0, paux.getBlocks()[-1].getSize()):
        dir.append(pts.Directions.U)
    """print(fold1, len(fold1))
    print(fold2, len(fold2))
    print(paux.getBlocks()[0]._seq)
    print(paux.getBlocks()[-1]._seq)
    print(sep.getSize())
    print(len(p), len(dir))"""
    return dir

def prot_coord(dir: list[pts.Directions], in3D = False):
    coord_x = [0]
    coord_y = [0]
    coord_z = [0]
    if in3D:
        ind_to_coord = {0 : [0,0,0]}
        coord_to_ind = {(0,0,0) : 0}
    else:
        ind_to_coord = {0 : [0,0]}
        coord_to_ind = {(0,0) : 0}
    i = 1
    for d in dir:
        if d == pts.Directions.U:
            coord_x.append(coord_x[-1])
            coord_y.append(coord_y[-1] + 1)
            coord_z.append(coord_z[-1])
        elif d == pts.Directions.D:
            coord_x.append(coord_x[-1])
            coord_y.append(coord_y[-1] - 1)
            coord_z.append(coord_z[-1])
        elif d == pts.Directions.L:
            coord_x.append(coord_x[-1] - 1)
            coord_y.append(coord_y[-1])
            coord_z.append(coord_z[-1])
        elif d == pts.Directions.R:
            coord_x.append(coord_x[-1] + 1)
            coord_y.append(coord_y[-1])
            coord_z.append(coord_z[-1])
        elif d == pts.Directions.F:
            coord_x.append(coord_x[-1])
            coord_y.append(coord_y[-1])
            coord_z.append(coord_z[-1] + 1)
        elif d == pts.Directions.B:
            coord_x.append(coord_x[-1])
            coord_y.append(coord_y[-1])
            coord_z.append(coord_z[-1] - 1)
        if in3D:
            ind_to_coord[i] = (coord_x[-1], coord_y[-1], coord_z[-1])
            if coord_to_ind.get((coord_x[-1], coord_y[-1], coord_z[-1])) != None:
                print('HAY ELEMENTOS COLINDANTES', (coord_x[-1], coord_y[-1], coord_z[-1]))
            coord_to_ind[(coord_x[-1], coord_y[-1], coord_z[-1])] = i
        else:
            ind_to_coord[i] = (coord_x[-1], coord_y[-1])
            if coord_to_ind.get((coord_x[-1], coord_y[-1])) != None:
                print('HAY ELEMENTOS COLINDANTES', (coord_x[-1], coord_y[-1]))
            coord_to_ind[(coord_x[-1], coord_y[-1])] = i
        i += 1
    coord = []
    if in3D:
        coord = [(coord_x[i], coord_y[i], coord_z[i]) for i in range(0, len(coord_x))]
    else:
        coord = [(coord_x[i], coord_y[i]) for i in range(0, len(coord_x))] 
    return ind_to_coord, coord_to_ind, coord

def next_dir(d: pts.Directions, c: tuple[int, int, int], in3D = True):
    if in3D:
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
    else:
        if d == pts.Directions.U:
            return (c[0], c[1] + 1)
        elif d == pts.Directions.D:
            return (c[0], c[1] - 1)
        elif d == pts.Directions.L:
            return (c[0] - 1, c[1])
        elif d == pts.Directions.R:
            return (c[0] + 1, c[1])

def neighbours(c: tuple[int, int, int], in3D = True) -> list[pts.Directions]:
    res = []
    if in3D:
        for d in pts.Directions:
            res.append(next_dir(d, c, in3D))
    else:
        for d in [dir for dir in pts.Directions if dir != pts.Directions.F and dir != pts.Directions.B]:
            res.append(next_dir(d, c, in3D))
    return res

def fitness(p: str, ind_to_coor: dict, coor_to_ind: dict, in3D = True) -> int:
    energy = 0
    for i in range(len(p)):
        for n in neighbours(ind_to_coor[i], in3D):
            if coor_to_ind.get(n) != None and coor_to_ind.get(n) > i and coor_to_ind.get(n) - i > 1 and p[i] == 'H' and p[coor_to_ind.get(n)] == 'H':
                energy += 4
            elif coor_to_ind.get(n) != None and coor_to_ind.get(n) > i and coor_to_ind.get(n) - i > 1 and (p[i] == 'P' and p[coor_to_ind.get(n)] == 'N' or p[i] == 'N' and p[coor_to_ind.get(n)] == 'P'):
                energy += 1
            elif coor_to_ind.get(n) != None and coor_to_ind.get(n) > i and coor_to_ind.get(n) - i > 1 and (p[i] == 'P' and p[coor_to_ind.get(n)] == 'P' or p[i] == 'N' and p[coor_to_ind.get(n)] == 'N'):
                energy -= 1
    return energy
        
def prot_fold(str_seq: str):
    res = aproxHPNX(str_seq)
    ind_to_coord, coord_to_ind, coord = prot_coord(res, True)
    fit = fitness(str_seq, ind_to_coord, coord_to_ind)
    print('El valor obtenido por el algoritmo aproximado es', fit)
    """color = []
    for ch in str_seq:
        if ch == 'H':
            color.append('black')
        elif ch == 'X':
            color.append('white')
        elif ch == 'P':
            color.append('blue')
        else:
            color.append('yellow')
    coord_x = [t[0] for t in coord]
    coord_y = [t[1] for t in coord] 
    coord_z = [t[2] for t in coord]
    ax = plt.axes(projection ='3d')
    ax.plot3D(coord_x, coord_y, coord_z, c = 'black', zorder = 1)
    ax.scatter3D(coord_x, coord_y, coord_z, c = color, zorder = 2)
    plt.axis('equal')
    plt.show()"""
    return coord, fit


"""for i in range(100):
    r = random.randint(200, 500)
    str_seq = ''
    for _ in range(r):
        str_seq += random.choice(['H', 'P', 'N', 'X'])
    print("Caso %i: %s" % (i, str_seq))
    prot_fold(str_seq)"""