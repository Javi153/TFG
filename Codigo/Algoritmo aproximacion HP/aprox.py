import proteins as pts
from matplotlib import pyplot as plt

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
    next = '0'
    cant = 0
    res = []
    for ch in str_seq:
        if ch == next:
            cant += 1
        else:
            res.append(cant)
            if next == '0':
                next = '1'
            else:
                next = '0'
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
            coord_to_ind[(coord_x[-1], coord_y[-1], coord_z[-1])] = i
        else:
            ind_to_coord[i] = (coord_x[-1], coord_y[-1])
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
    paux = pts.prot(p, False)
    pb = pts.prot_block(pts.Block_type.SEP)
    z_i = 0
    pb.add_amino(0)
    for _ in range(0, paux.getBlocks()[0].getSize()):
        dir.append(pts.Directions.D)
    P1, P2, reverse = subroutine1(pts.prot([pb] + paux.getBlocks()[1:-1] + [pb]))
    if reverse:
        z_i = P2.getBlocks()[-1].getSize()
    else:
        z_i = P1.getBlocks()[-1].getSize()
    if P1.Ny() > P2.Nx():
        if reverse:
            z_i += P1.del_first_h()
        else:
            z_i += P1.del_last_h()
    elif P1.Ny() < P2.Nx():
        if reverse:
            z_i += P2.del_last_h()
        else:
            z_i += P2.del_first_h()
    for _ in range(z_i // 2):
        dir.append(pts.Directions.D)
    dir.append(pts.Directions.R)
    for _ in range(z_i // 2):
        dir.append(pts.Directions.U)
    fold1 = P1.fold(pts.Block_type.X_BLOCK, reverse)
    fold2 = P2.fold(pts.Block_type.Y_BLOCK, not reverse)
    if P1.Ny() == P2.Nx():
        if reverse:
            fold2[0] = pts.Directions.L
        else:
            fold1[0] = pts.Directions.L
    for _ in range(0, paux.getBlocks()[-1].getSize()):
        dir.append(pts.Directions.U)
    return dir

def algorithmC(p: list[int]) -> list[pts.Directions]:
    dir = []
    return dir

def prot_fold(str_seq: str, algorithm: str):
    seq = format_seq(str_seq)
    res = []
    if algorithm == 'A':
        res = algorithmA(seq)
    elif algorithm == 'B':
        res = algorithmB(seq)
    else:
        res = algorithmC(seq)
    if algorithm == 'A' or algorithm == 'B':
        _, _, coord = prot_coord(res)
        color = []
        for ch in str_seq:
            if ch == '1':
                color.append('black')
            else:
                color.append('white')
        coord_x = [t[0] for t in coord]
        coord_y = [t[1] for t in coord] 
        plt.plot(coord_x, coord_y, '-', c = 'black', zorder = 1)
        plt.axis('equal')
        plt.scatter(coord_x, coord_y, c = color, s = 100, edgecolors = 'black', zorder = 2)
        plt.grid(color = 'black', linewidth=0.5)
        plt.show()
    return coord

    
str_seq = '0101000111100100000001101011000101101001'
prot_fold(str_seq, 'A')

str_seq = '0101101001001011101001101001010'
prot_fold(str_seq, 'A')

str_seq = '00100010001000100100010001000100'
prot_fold(str_seq, 'A')

str_seq = '1'
prot_fold(str_seq, 'A')