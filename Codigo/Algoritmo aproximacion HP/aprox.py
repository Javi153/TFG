import proteins as pts
from matplotlib import pyplot as plt
import math
from mpl_toolkits import mplot3d
import random
import pandas as pd
import numpy as np

def Mxy(p1: pts.prot, p2: pts.prot, getmin = True) -> int:
    if getmin:
        return min(p1.Nx(), p2.Ny())
    else:
        return max(p1.Nx(), p2.Ny())
    
def Myx(p1: pts.prot, p2: pts.prot, getmin = True) -> int:
    return Mxy(p2, p1, getmin)
    
def cond(a, b, A, B) -> bool:
    return a < b or (a == b and A < B)

def subroutine1(p : pts.prot):
    if len(p.getBlocks()) < 3: 
        return p, pts.prot([], False), False
    reverse = False
    P1 = []
    P2 = []
    pb = pts.prot_block(pts.Block_type.SEP)
    pb.add_amino(0)
    B1 = pts.prot(p.getBlocks()[0:2] + [pb])
    sep = p.getBlocks()[2]
    B2 = pts.prot([pb] + p.getBlocks()[3:])
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
        B1 = pts.prot(p.getBlocks()[0:2*(i-1)] + [pb])
        B2 = pts.prot([pb] + p.getBlocks()[(2*i-1):])
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
    next = 'P'
    cant = 0
    res = []
    for ch in str_seq:
        if ch == next:
            cant += 1
        else:
            res.append(cant)
            if next == 'P':
                next = 'H'
            else:
                next = 'P'
            cant = 1
    if cant > 0:
        res.append(cant)
    if len(res) % 2 == 0:
        res.append(0)
    return res

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
                print('HAY ELEMENTOS COLINDANTES')
            coord_to_ind[(coord_x[-1], coord_y[-1], coord_z[-1])] = i
        else:
            ind_to_coord[i] = (coord_x[-1], coord_y[-1])
            if coord_to_ind.get((coord_x[-1], coord_y[-1])) != None:
                print('HAY ELEMENTOS COLINDANTES')
            coord_to_ind[(coord_x[-1], coord_y[-1])] = i
        i += 1
    coord = []
    if in3D:
        coord = [(coord_x[i], coord_y[i], coord_z[i]) for i in range(0, len(coord_x))]
    else:
        coord = [(coord_x[i], coord_y[i]) for i in range(0, len(coord_x))] 
    return ind_to_coord, coord_to_ind, coord


def algorithmA(p: list[int]) -> list[pts.Directions]:
    dir = []
    if sum(p) <= 3 or len(p) == 1:
        for _ in range(sum(p) - 1):
            dir.append(pts.Directions.R)
        return dir
    paux = pts.prot(p, False)
    pb = pts.prot_block(pts.Block_type.SEP)
    pb.add_amino(0)
    for _ in range(0, paux.getBlocks()[0].getSize()):
        dir.append(pts.Directions.D)
    P1, sep, P2, reverse = subroutine1(pts.prot([pb] + paux.getBlocks()[1:-1] + [pb]))
    fold1 = P1.fold(pts.Block_type.Y_BLOCK, reverse)
    fold2 = P2.fold(pts.Block_type.X_BLOCK, not reverse)
    union = []
    for _ in range(sep.getSize() // 2):
        union.append(pts.Directions.D)
    if fold1 and fold2:
        union.append(pts.Directions.R)
    for _ in range(sep.getSize() // 2):
        union.append(pts.Directions.U)
    if reverse:
        dir = dir + fold2 + union + fold1
    else:
        dir = dir + fold1 + union + fold2
    for _ in range(0, paux.getBlocks()[-1].getSize()):
        dir.append(pts.Directions.U)
    return dir

def algorithmB(p: list[int]) -> list[pts.Directions]:
    dir = []
    if sum(p) <= 3 or len(p) == 1:
        for _ in range(sum(p) - 1):
            dir.append(pts.Directions.R)
        return dir
    paux = pts.prot(p, False)
    pb = pts.prot_block(pts.Block_type.SEP)
    z_i = 0
    pb.add_amino(0)
    P1, sep, P2, reverse = subroutine1(pts.prot([pb] + paux.getBlocks()[1:-1] + [pb]))
    if len(P1.getBlocks()) <= 1 or len(P2.getBlocks()) <= 1:
        return [pts.Directions.R for _ in range(sum(p) - 1)]
    z_i = sep.getSize()
    if P1.Ny() > P2.Nx() and z_i == 0:
        if reverse:
            z_i += P1.del_first_from_block(pts.Block_type.Y_BLOCK)
        else:
            z_i += P1.del_last_from_block(pts.Block_type.Y_BLOCK)
    elif P1.Ny() < P2.Nx() and z_i == 0:
        if reverse:
            z_i += P2.del_last_from_block(pts.Block_type.X_BLOCK)
        else:
            z_i += P2.del_first_from_block(pts.Block_type.X_BLOCK)
    fold1 = P1.fold(pts.Block_type.Y_BLOCK, reverse)
    fold2 = P2.fold(pts.Block_type.X_BLOCK, not reverse)
    union = []
    for _ in range(z_i // 2):
        union.append(pts.Directions.D)
    if z_i > 1 or (fold1 or fold2):
        union.append(pts.Directions.R)
    for _ in range(z_i // 2):
        union.append(pts.Directions.U)
    if reverse:
        dir = fold2 + union + fold1
    else:
        dir = fold1 + union + fold2
    if not fold1 or not fold2:
        if dir == []:
            dir = [pts.Directions.R]
        return [pts.Directions.D for _ in range(paux.getBlocks()[0].getSize())] + dir + [pts.Directions.U for _ in range(paux.getBlocks()[-1].getSize())]
    if P1.Ny() == P2.Nx():
        if P1.Ny() != 1:
            first_one = 0
            pre_last_one = 0
            post_last_one = 0
            first_prot = P1
            last_prot = P2
            first_type = pts.Block_type.Y_BLOCK
            last_type = pts.Block_type.X_BLOCK
            if reverse:
                first_prot, last_prot, first_type, last_type = last_prot, first_prot, last_type, first_type
            for bl in first_prot.getBlocks():
                if bl.getType() != first_type:
                    first_one += bl.getSize()
                else:
                    break
            for i in range(len(last_prot.getBlocks()) - 1, -1, -1):
                if last_prot.getBlocks()[i].getType() != last_type:
                    pre_last_one += last_prot.getBlocks()[i].getSize()
                else:
                    post_last_one = last_prot.del_last_from_block(last_type) - pre_last_one
                    break
            for _ in range(first_one + 1):
                dir.pop(0)
            dir.insert(0, pts.Directions.L)
            for _ in range(first_one):
                dir.insert(0, pts.Directions.D)
            for _ in range(pre_last_one + post_last_one):
                dir.pop()
            dir.append(pts.Directions.R)
            for _ in range(post_last_one // 2 - 1):
                dir.append(pts.Directions.R)
            if post_last_one != 1:
                dir.append(pts.Directions.U)
            for _ in range(post_last_one // 2 - 1):
                dir.append(pts.Directions.L)
            for _ in range(pre_last_one):
                dir.append(pts.Directions.U)
        dir = [pts.Directions.D for _ in range(paux.getBlocks()[0].getSize())] + dir + [pts.Directions.U for _ in range(paux.getBlocks()[-1].getSize())]
    else:
        first = P1
        second = P2
        shorter = P1
        first_type = pts.Block_type.Y_BLOCK
        second_type = pts.Block_type.X_BLOCK
        if reverse:
            first, second, first_type, second_type = second, first, second_type, first_type
        if P2.Nx() < P1.Ny():
            shorter = P2
        if shorter == first:
            hor1 = pts.Directions.R
            hor2 = pts.Directions.L
        else:
            hor1 = pts.Directions.L
            hor2 = pts.Directions.R
        last_one = 0
        second_last_one = 0
        aux = 0
        while abs(P1.Ny() - P2.Nx()) != 1:
            if shorter == first:
                aux = second.del_last_from_block(second_type)
            else:
                aux = first.del_first_from_block(first_type)
            second_last_one = last_one
            last_one += aux
        if last_one == 0:
            if shorter == first:
                for bl in second.getBlocks()[::-1]:
                    if bl.getType() != second_type:
                        last_one += bl.getSize()
                    else:
                        break
            else:
                for bl in first.getBlocks():
                    if bl.getType() != first_type:
                        last_one += bl.getSize()
                    else:
                        break
        if shorter != first:
            for_counter = range(last_one, -1, -1)
        else:
            for_counter = range(len(dir) - last_one - 1, len(dir))
        for i in for_counter:
            if dir[i] == pts.Directions.U:
                dir[i] = hor2
            elif dir[i] == hor2:
                dir[i] = pts.Directions.D
            elif dir[i] == pts.Directions.D:
                dir[i] = hor1
            elif dir[i] == hor1:
                dir[i] = pts.Directions.U
        cond = (P1.Ny() - P2.Nx() > 1 and ((reverse and p[0] == 0 and P2.getBlocks()[1].getType() == pts.Block_type.X_BLOCK) or (not reverse and p[-1] == 0 and P2.getBlocks()[-2].getType() == pts.Block_type.X_BLOCK))) or (P2.Nx() - P1.Ny() > 1 and ((reverse and p[-1] == 0 and P1.getBlocks()[-2].getType() == pts.Block_type.Y_BLOCK) or (not reverse and p[0] == 0 and P1.getBlocks()[1].getType() == pts.Block_type.Y_BLOCK)))
        if cond:
            cont = last_one - second_last_one
            if shorter == first:
                for i in range(cont // 2 - 1):
                    dir[len(dir) - last_one + i] = pts.Directions.U
                dir[len(dir) - last_one + cont // 2 - 1] = pts.Directions.L
                for i in range(cont // 2):
                    dir[len(dir) - last_one + cont // 2 + i] = pts.Directions.D
                dir[len(dir) - 1 - second_last_one] = pts.Directions.L
                """dir[len(dir) - second_last_one - 1] = pts.Directions.D
                if second_last_one != 0 and dir[len(dir) - second_last_one] == pts.Directions.U:
                    dir[len(dir) - second_last_one] = pts.Directions.L
                    for i in range(len(dir) - second_last_one + 1, len(dir)):
                        if dir[i] == pts.Directions.D:
                            dir[i] = pts.Directions.L
                            break"""
            else:
                for i in range(cont // 2):
                    dir[second_last_one + 1 + i] = pts.Directions.U
                dir[second_last_one + cont // 2 + 1] = pts.Directions.L
                for i in range(cont // 2 - 1):
                    dir[second_last_one + cont // 2 + 2 + i] = pts.Directions.D
                dir[second_last_one] = pts.Directions.L
                """dir[second_last_one] = pts.Directions.U
                if second_last_one != 0 and dir[second_last_one - 1] == pts.Directions.D:
                    dir[second_last_one - 1] = pts.Directions.L
                    for i in range(second_last_one - 2, -1, -1):
                        if dir[i] == pts.Directions.U:
                            dir[i] = pts.Directions.L
                            break"""
        else:
            aux = 0
            if shorter == first:
                for bl in first.getBlocks():
                    if bl.getType() != first_type:
                        aux += bl.getSize()
                    else:
                        break
                for i in range(0, aux):
                    dir[i] = pts.Directions.R
            else:
                for bl in second.getBlocks()[::-1]:
                    if bl.getType() != second_type:
                        aux += bl.getSize()
                    else:
                        break
                for i in range(len(dir) - 1, len(dir) - 1 - aux, -1):
                    dir[i] = pts.Directions.R
        dir = [hor1 for _ in range(paux.getBlocks()[0].getSize())] + dir + [hor2 for _ in range(paux.getBlocks()[-1].getSize())]    
        """
        if shorter != first:
            while P1.Ny() != P2.Nx():
                aux = first.del_first_h()
                pre_last_one += aux
                post_last_one = aux
            pre_last_one -= post_last_one
            for _ in range(pre_last_one + post_last_one):
                dir.pop(0)
            for _ in range((post_last_one // 2) - 1):
                dir.insert(0, pts.Directions.R)
            if post_last_one == 1:
                dir.insert(0, pts.Directions.L)
            else:
                dir.insert(0, pts.Directions.D)
            for _ in range(post_last_one // 2):
                dir.insert(0, pts.Directions.L)
            for _ in range(pre_last_one):
                dir.insert(0, pts.Directions.L)
        else:
            while P1.Ny() != P2.Nx():
                aux = second.del_last_h()
                pre_last_one += aux
                post_last_one = aux
            pre_last_one -= post_last_one
            for _ in range(pre_last_one + post_last_one):
                dir.pop()
            for _ in range((post_last_one // 2) - 1):
                dir.append(pts.Directions.R)
            if post_last_one == 1:
                dir.append(pts.Directions.L)
            else:
                dir.append(pts.Directions.U)
            for _ in range(post_last_one // 2):
                dir.append(pts.Directions.L)
            for _ in range(pre_last_one):
                dir.append(pts.Directions.L)
        dir = [hor1 for _ in range(paux.getBlocks()[0].getSize())] + dir + [hor2 for _ in range(paux.getBlocks()[-1].getSize())]    
        """
    return dir

def algorithmC(p: list[int], f) -> list[pts.Directions]:
    dir = []
    if sum(p) <= 3 or len(p) == 1:
        for _ in range(sum(p) - 1):
            dir.append(pts.Directions.R)
        return dir
    paux = pts.prot(p, False)
    pb = pts.prot_block(pts.Block_type.SEP)
    pb.add_amino(0)
    P1, sep, P2, reverse = subroutine1(pts.prot([pb] + paux.getBlocks()[1:-1] + [pb]))
    K = math.floor(f(min(P1.Ny(), P2.Nx())))
    J = math.floor((min(P1.Ny(), P2.Nx()) - 2 * K + 1) / K)
    if K <= 1 or J < 2:
        res = algorithmB(p)
        return res
    dir.append(pts.Directions.R)
    for _ in range(sep.getSize() // 2):
        dir.insert(0, pts.Directions.D)
        dir.append(pts.Directions.U)
    first_p = P1
    second_p = P2
    first = P1.fold(pts.Block_type.Y_BLOCK, reverse)
    second = P2.fold(pts.Block_type.X_BLOCK, not reverse)
    first_type = pts.Block_type.Y_BLOCK
    second_type = pts.Block_type.X_BLOCK
    if reverse:
        first_p, second_p = second_p, first_p
        first, second = second, first
        first_type, second_type = second_type, first_type
    ind = len(first) - 1
    h = 0
    for t in range(1, K + 1):
        while h < 2 * (J - 1):
            if first[ind] == pts.Directions.D:
                h += 1
                if h % 2 == 0:
                    first_p.del_last_from_block(first_type)
                if t % 2 == 0:
                    first[ind] = pts.Directions.U
            elif first[ind] == pts.Directions.R and t % 2 == 0:
                first[ind] = pts.Directions.L
            elif first[ind] == pts.Directions.L and t % 2 == 0:
                first[ind] = pts.Directions.R
            ind -= 1
        if h == 2 * (J - 1) and t != K:
            if t % 2 != 0:
                s = first_p.del_last_from_block(first_type)
                for _ in range(s // 2 - 1):
                    first[ind] = pts.Directions.R
                    ind -= 1
                first[ind] = pts.Directions.D
                ind -= 1
                for _ in range(s // 2):
                    first[ind] = pts.Directions.L
                    ind -= 1
                s = first_p.del_last_from_block(first_type)
                for _ in range(s // 2 - 1):
                    first[ind] = pts.Directions.D
                    ind -= 1
                first[ind] = pts.Directions.L
                ind -= 1
                for _ in range(s // 2):
                    first[ind] = pts.Directions.U
                    ind -= 1
                first[ind] = pts.Directions.F
                ind -= 1
                s = 0
                while(s <= 2):
                    s += first_p.del_last_from_block(first_type)
                for _ in range(s // 2 - 2):
                    first[ind] = pts.Directions.L
                    ind -= 1
                first[ind] = pts.Directions.U
                ind -= 1
                for _ in range(s // 2 - 1):
                    first[ind] = pts.Directions.R
                    ind -= 1
                first[ind] = pts.Directions.U
                ind -= 1
            else:
                s = 0
                while s <= 2:
                    s += first_p.del_last_from_block(first_type)
                first[ind] = pts.Directions.F
                ind -= 1
                for _ in range(s // 2 - 1):
                    first[ind] = pts.Directions.R
                    ind -= 1
                first[ind] = pts.Directions.D
                ind -= 1
                for _ in range(s // 2 - 2):
                    first[ind] = pts.Directions.L
                    ind -= 1
                first[ind] = pts.Directions.D
                ind -= 1
        h = 2
    if K % 2 == 0:
        while ind >= 0:
            if first[ind] == pts.Directions.D:
                first[ind] = pts.Directions.U
            elif first[ind] == pts.Directions.R:
                first[ind] = pts.Directions.L
            elif first[ind] == pts.Directions.L:
                first[ind] = pts.Directions.R
            ind -= 1
    ind = 0
    h = 0
    for t in range(1, K + 1):
        while h < 2 * (J - 1):
            if second[ind] == pts.Directions.U:
                h += 1
                if h % 2 == 0:
                    second_p.del_first_from_block(second_type)
                if t % 2 == 0:
                    second[ind] = pts.Directions.D
            elif second[ind] == pts.Directions.R and t % 2 == 0:
                second[ind] = pts.Directions.L
            elif second[ind] == pts.Directions.L and t % 2 == 0:
                second[ind] = pts.Directions.R
            ind += 1
        if h == 2 * (J - 1) and t != K:
            if t % 2 == 0:
                s = second_p.del_first_from_block(second_type)
                for _ in range(s // 2 - 1):
                    second[ind] = pts.Directions.L
                    ind += 1
                second[ind] = pts.Directions.D
                ind += 1
                for _ in range(s // 2):
                    second[ind] = pts.Directions.R
                    ind += 1
                s = second_p.del_first_from_block(second_type)
                for _ in range(s // 2 - 1):
                    second[ind] = pts.Directions.D
                    ind += 1
                second[ind] = pts.Directions.R
                ind += 1
                for _ in range(s // 2):
                    second[ind] = pts.Directions.U
                    ind += 1
                second[ind] = pts.Directions.B
                ind += 1
                s = 0
                while s <= 2:
                    s += second_p.del_first_from_block(second_type)
                for _ in range(s // 2 - 2):
                    second[ind] = pts.Directions.R
                    ind += 1
                second[ind] = pts.Directions.U
                ind += 1
                for _ in range(s // 2 - 1):
                    second[ind] = pts.Directions.L
                    ind += 1
                second[ind] = pts.Directions.U
                ind += 1
            else:
                s = 0
                while s <= 2:
                    s += second_p.del_first_from_block(second_type)
                second[ind] = pts.Directions.B
                ind += 1
                for _ in range(s // 2 - 1):
                    second[ind] = pts.Directions.L
                    ind += 1
                second[ind] = pts.Directions.D
                ind += 1
                for _ in range(s // 2 - 2):
                    second[ind] = pts.Directions.R
                    ind += 1
                second[ind] = pts.Directions.D
                ind += 1
        h = 2
    if K % 2 == 0:
        while ind < len(second):
            if second[ind] == pts.Directions.U:
                second[ind] = pts.Directions.D
            elif second[ind] == pts.Directions.R:
                second[ind] = pts.Directions.L
            elif second[ind] == pts.Directions.L:
                second[ind] = pts.Directions.R
            ind += 1
    if K % 2 == 0:
        dir = [pts.Directions.U for _ in range(paux.getBlocks()[0].getSize())] + first + dir + second + [pts.Directions.D for _ in range(paux.getBlocks()[-1].getSize())]    
    else:
        dir = [pts.Directions.D for _ in range(paux.getBlocks()[0].getSize())] + first + dir + second + [pts.Directions.U for _ in range(paux.getBlocks()[-1].getSize())]            
    return dir

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
                energy += 1
    return energy

def prot_fold(str_seq: str, algorithm: str, f = None):
    print(str_seq)
    seq = format_seq(str_seq)
    res = []
    if algorithm == 'A':
        res = algorithmA(seq)
    elif algorithm == 'B':
        res = algorithmB(seq)
    elif algorithm == 'C' and f != None:
        res = algorithmC(seq, f)
    else:
        print('El algoritmo especificado no existe o no se pasó una función al algoritmo C')
        return
    color = []
    for ch in str_seq:
        if ch == 'H':
            color.append('black')
        else:
            color.append('white')
    if algorithm == 'A' or algorithm == 'B':
        ind_to_coor, coor_to_ind, coord = prot_coord(res)
        print('El valor obtenido por el algoritmo %s es' % algorithm, fitness(str_seq, ind_to_coor, coor_to_ind, False))
        coord_x = [t[0] for t in coord]
        coord_y = [t[1] for t in coord] 
        plt.plot(coord_x, coord_y, '-', c = 'black', zorder = 1)
        plt.axis('equal')
        plt.scatter(coord_x, coord_y, c = color, s = 100, edgecolors = 'black', zorder = 2)
        #plt.grid(color = 'gray', linewidth=0.5)
        plt.show()
    else:
        ind_to_coor, coor_to_ind, coord = prot_coord(res, in3D = True)
        print('El valor obtenido por el algoritmo C es', fitness(str_seq, ind_to_coor, coor_to_ind))
        coord_x = [t[0] for t in coord]
        coord_y = [t[1] for t in coord] 
        coord_z = [t[2] for t in coord]
        ax = plt.axes(projection ='3d')
        ax.plot3D(coord_x, coord_y, coord_z, c = 'black', zorder = 1)
        ax.scatter3D(coord_x, coord_y, coord_z, c = color, zorder = 2, s = 30)
        plt.axis('equal')
        plt.show()
    return coord

"""def trans_amino(amino: str) -> str:
    if amino in {"GLY", "ALA", "VAL", "LEU", "ILE", "MET", "PHE", "TRP", "PRO", "G", "A", "V", "L", "I", "M", "F", "W", "P"}:
        return 'H'
    else:
        return 'P'
    
data_file = open("prueba.txt",'r')
for line in data_file:
    data = line.split()

str_seq = ''
for item in data:
    str_seq += trans_amino(item)

aux = prot_fold(str_seq, 'C', math.sqrt)

result = open("result.txt", 'w')
result.write("MOLECULE algo\n")
i = 1
for coord in aux:
    result.write("ATOM %6i  CA %s A %4i     %6.3f  %6.3f  %6.3f  1.00  0.00  C\n" % (i, data[i-1], i, coord[0], coord[1], coord[2]))
    i += 1
result.write("END")
result.close()
data_file.close()"""

#str_seq = 'PHPHPPPHHHHPPHPPPPPPPHHPHPHHPPPHPHHPHPPH'
#prot_fold(str_seq, 'A')

#str_seq = 'PHPHHPHPPHPPHPHHHPHPPHHPHPPHPHP'
#prot_fold(str_seq, 'A')

#str_seq = 'PPHPPPHPPPHPPPHPPPHHPPPHPPPHPPPHPP'
#prot_fold(str_seq, 'B')

#str_seq = 'PHPHPPPHPPPHPPPHPPPHHPPPHPPPHPPPHPHP'
#prot_fold(str_seq, 'B')

#str_seq = 'PHPPPHPPPHPPPHPPHPPPHPPPHPPPHPPPHP'
#prot_fold(str_seq, 'B')

#str_seq = 'PHPPHPPHPPPHPPPHPPPHPPHPPPHPPPHPPPHP'
#prot_fold(str_seq, 'B')

#str_seq = 'PHPHPPPHHHHPPHPPPPPPPHHPHPHHPPPHPHHPHPPH'
#prot_fold(str_seq, 'B')

#str_seq = 'HPPPHPPPHPPPHHPPPHPPPHPPPHPPPHPPPHPPPHP'
#prot_fold(str_seq, 'A')

#str_seq = 'HPPPHPPPHPPPHPPPPHPPPHPPPHPPPHPPPHPPHPP'
#prot_fold(str_seq, 'B')

str_seq = 'HPPPHPPPHPPPHHPPPHPPPHPPPHPPPHPPPHPPPHPHPPPHPPPHPPPHHPPPHPPPHPPPHPPPHPPPHPPPHPHPPPHPPPHPPPHHPPPHPPPHPPPHPPPHPPPHPPPHPHPPPHPPPHPPPHHPPPHPPPHPPPHPPPHPPPHPPPHP'
prot_fold(str_seq, 'C', math.sqrt)

#str_seq = 'HPPHPHHPHHPPHPHPHHHPHPHPPHPHPHHPHPPHHHPHPHHPHPPPPHHHPPHHPHPHHPHHPPHPHPPHPPHHHPPHPPHPPHPHHPHPHPPHPHPHHPHPPHHHPHPHHPHPPPPHHHPPHHPHPHHPHH'
#prot_fold(str_seq, 'C', math.sqrt)

#str_seq = 'HPPHHHPHPPPHHPHPPPHPHHHPPHHPHPHPHPHHHPPHPHPPHPHPHHHPHPHPHHHPHPPPHHPPHPPPHPPHPHPPPPHHPHPPHPHPHHPPHPHPPHPHHPPHHPHPHHPHPHPPHHHPHP'
#prot_fold(str_seq, 'C', math.sqrt)

#str_seq = 'PPHPPPHPPPHPPPHPPHPPPPPPPHPHHHPP'
#prot_fold(str_seq, 'A')

#str_seq = 'PPPPPHHPHPPPHPPPPPHHHPPHHHHPPPHHPPHPPHPHHHHPHHPHPPHPHPHHHPHHPH'
#prot_fold(str_seq, 'C', math.sqrt)

#str_seq = 'PHHPPHHPHHHPHPPHHPPHPPPPPHHHHHHHHHPHHPPHHHPPPPP'
#prot_fold(str_seq, 'C', math.sqrt)

"""for i in range(100):
    r = random.randint(1, 50)
    str_seq = ''
    for _ in range(r):
        str_seq += random.choice(['H', 'P'])
    print("Caso %i: %s" % (i, str_seq))
    prot_fold(str_seq, 'B')"""