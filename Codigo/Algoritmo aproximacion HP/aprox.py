from enum import Enum
from matplotlib import pyplot as plt
from cycler import cycler

Block_type = Enum('Block_type', ['SEP', 'X_BLOCK', 'Y_BLOCK'])

class prot_block:
    _type = Block_type.SEP
    _seq = []
    _num_aminos = 0
    _num_h = 0
    
    def __init__(self, type, seq):
        self._type = type
        self._seq = seq
        self._num_aminos = sum(seq)
        if type == Block_type.SEP:
            self._num_h = 0
        else:
            self._num_h = sum([seq[i] for i in range(0, len(seq)) if i % 2 == 0])
        
    def __init__(self, type):
        self._seq = []
        self._type = type
        self._num_aminos = 0
        self._num_h = 0
        
    def add_amino(self, am):
        self._seq.append(am)
        self._num_aminos += am
        if len(self._seq) % 2 != 0 and self._type != Block_type.SEP:
            self._num_h += am
    
    def del_amino(self):
        if self._seq:
            aux = self._seq.pop()
            self._num_aminos -= aux
            if len(self._seq) %2 == 0 and self._type != Block_type.SEP:
                self._num_h -= aux
                
    def del_last_h(self):
        if self._seq:
            aux = self._seq.pop()
            self._num_aminos -= aux
            if len(self._seq) %2 != 0:
                aux = self._seq.pop()
                self._num_aminos -= aux
            self._num_h -= aux
    
    def del_first_h(self):
        if self._seq:
            aux = self._seq.pop()
            self._num_aminos -= aux
            if len(self._seq) %2 != 0:
                aux = self._seq.pop()
                self._num_aminos -= aux
            self._num_h -= aux

    def fold(self, ftype, reverse = False):
        if self._type == ftype:
            if reverse:
                vert = Directions.U
                hor = Directions.R
                hor2 = Directions.L
            else:
                vert = Directions.D
                hor = Directions.L
                hor2 = Directions.R
            dir = []
            for i in range(0, len(self._seq)):
                if i % 2 == 0 and self._seq[i] == 1:
                    dir.append(vert)
                elif self._seq[i] != 0:
                    for _ in range(0, self._seq[i] // 2):
                        dir.append(hor)
                    dir.append(vert)
                    for _ in range(0, self._seq[i] // 2):
                        dir.append(hor2)
            return dir
        else:
            return []
                
    def getType(self):
        return self._type
    
    def getN(self):
        return self._num_h
    
    def getSize(self):
        return self._num_aminos
    
    def setType(self, type1):
        self._type = type1
    
class prot:
    _blocks : list[prot_block]
    _nx = 0
    _ny = 0
        
    def __init__(self, blocks, isBlock = True):
        self._nx = 0
        self._ny = 0
        if isBlock:
            self._blocks = blocks
            for bl in blocks:
                if bl.getType() == Block_type.X_BLOCK:
                    self._nx += bl.getN()
                elif bl.getType() == Block_type.Y_BLOCK:
                    self._ny += bl.getN()
        else:
            self.update_partition(blocks)

    def del_amino(self):
        if self._blocks[len(self._blocks) - 1].getN() == 1:
            self._blocks.pop()
        else:
            self._blocks[len(self._blocks) - 1].del_amino()
            
    def del_last_h(self):
        if self._blocks[len(self._blocks) - 1].getN() <= 1:
            self._blocks.pop()
        else:
            self._blocks[len(self._blocks) - 1].del_last_h()
            
    def del_first_h(self):
        if self._blocks[0].getN() <= 1:
            self._blocks.pop(0)
        else:
            self._blocks[0].del_first_h()

    def update_partition(self, seq):
        next_y = True
        carry = False
        bt = Block_type.Y_BLOCK
        pb = prot_block(Block_type.SEP)
        pb.add_amino(seq[0])
        self._blocks = []
        self._blocks.append(pb)
        i = 1
        while i < len(seq):
            if next_y:
                bt = Block_type.Y_BLOCK
            else:
                bt = Block_type.X_BLOCK
            next_y = not next_y
            pb = prot_block(bt)
            if carry:
                carry = False
                pb.add_amino(1)
                if bt == Block_type.X_BLOCK:
                    self._nx += 1
                else:
                    self._ny += 1
                i += 1
            while i < len(seq) and ((i % 2 != 0 and seq[i] == 1) or (i % 2 == 0 and seq[i] % 2 != 0)):
                pb.add_amino(seq[i])
                if i % 2 != 0:
                    if bt == Block_type.X_BLOCK:
                        self._nx += 1
                    else:
                        self._ny += 1
                i += 1
            if i % 2 == 0:
                self._blocks.append(pb)
                if i < len(seq):
                    pb = prot_block(Block_type.SEP)
                    pb.add_amino(seq[i])
                    self._blocks.append(pb)
                    i += 1
            elif i < len(seq):
                pb.add_amino(1)
                if bt == Block_type.X_BLOCK:
                    self._nx += 1
                else:
                    self._ny += 1
                self._blocks.append(pb)
                pb = prot_block(Block_type.SEP)
                pb.add_amino(0)
                self._blocks.append(pb)
                for _ in range(0, seq[i] - 2):
                    if next_y:
                        bt = Block_type.Y_BLOCK
                    else:
                        bt = Block_type.X_BLOCK
                    next_y = not next_y
                    pb = prot_block(bt)
                    pb.add_amino(1)
                    if bt == Block_type.X_BLOCK:
                        self._nx += 1
                    else:
                        self._ny += 1
                    self._blocks.append(pb)
                    pb = prot_block(Block_type.SEP)
                    pb.add_amino(0)
                    self._blocks.append(pb)
                carry = True
        if len(self._blocks) % 2 == 0:
            pb = prot_block(Block_type.SEP)
            pb.add_amino(0)
            self._blocks.append(pb)
        #REVISAR ESTA CONDICION
        if self._nx > self._ny or (self._nx == self._ny and self._blocks[1].getType() != Block_type.X_BLOCK and self._blocks[-2].getType() != Block_type.X_BLOCK):
            for bl in self._blocks:
                if bl.getType() == Block_type.X_BLOCK:
                    bl.setType(Block_type.Y_BLOCK)
                elif bl.getType() == Block_type.Y_BLOCK:
                    bl.setType(Block_type.X_BLOCK)
            self._nx, self._ny = self._ny, self._nx
            
    def Nx(self):
        return self._nx

    def Ny(self):
        return self._ny
    
    def getBlocks(self):
        return self._blocks
    
    def fold(self, ftype, reverse = False):
        if reverse:
            vert = Directions.U
            hor = Directions.R
            hor2 = Directions.L
        else:
            vert = Directions.D
            hor = Directions.L
            hor2 = Directions.R
        size = 0
        dir = []
        for bl in self._blocks:
            if bl.getType() == ftype:
                if size > 0:
                    for _ in range (0, size // 2):
                        dir.append(hor)
                    dir.append(vert)
                    for _ in range (0, size // 2):
                        dir.append(hor2)
                dir = dir + bl.fold(ftype, reverse)
                size = 0
            else:
                size += bl.getSize()
        if size > 0:
            for _ in range (0, size):
                dir.append(vert)
        return dir[:-1]
    
def Mxy(p1: prot, p2: prot, getmin = True) -> int:
    if getmin:
        return min(p1.Nx(), p2.Ny())
    else:
        return max(p1.Nx(), p2.Ny())
    
def Myx(p1: prot, p2: prot, getmin = True) -> int:
    return Mxy(p2, p1, getmin)
    
def cond(a, b, A, B) -> bool:
    return a < b or (a == b and A < B)

def subroutine1(p : prot):
    P1 = []
    P2 = []
    pb = prot_block(Block_type.SEP)
    pb.add_amino(0)
    if len(p.getBlocks()) > 1:
        B1 = prot(p.getBlocks()[0:3])
        B2 = prot([pb] + p.getBlocks()[3:])
        m1 = Mxy(B1, B2)
        m2 = Myx(B1, B2)
        M1 = Mxy(B1, B2, False)
        M2 = Myx(B1, B2, False)
        if m1 > m2:
            P1 = B2
            P2 = B1
            e = m1
            E = M1
        else:
            P1 = B1
            P2 = B2
            e = m2
            E = M2
        for i in range(3, ((len(p.getBlocks()) + 1) // 2) + 1):
            B1 = prot(p.getBlocks()[0:((i * 2) - 1)])
            B2 = prot([pb] + p.getBlocks()[((i * 2) - 1):])
            m1 = Mxy(B1, B2)
            m2 = Myx(B1, B2)
            M1 = Mxy(B1, B2, False)
            M2 = Myx(B1, B2, False)
            if cond(e, m1, E, M1):
                P1 = B2
                P2 = B1
                e = m1
                E = M1
            if cond(e, m2, E, M2):
                P1 = B1
                P2 = B2
                e = m2
                E = M2
    return P1, P2

Directions = Enum('Directions', ['U', 'D', 'L', 'R', 'F', 'B'])

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

def prot_coord2D(dir: list[Directions]):
    coord_x = [0]
    coord_y = [0]
    ind_to_coord = {0 : [0,0]}
    coord_to_ind = {(0,0) : 0}
    i = 1
    for d in dir:
        if d == Directions.U:
            coord_x.append(coord_x[len(coord_x)-1])
            coord_y.append(coord_y[len(coord_y)-1] + 1)
        elif d == Directions.D:
            coord_x.append(coord_x[len(coord_x)-1])
            coord_y.append(coord_y[len(coord_y)-1] - 1)
        elif d == Directions.L:
            coord_x.append(coord_x[len(coord_x)-1] - 1)
            coord_y.append(coord_y[len(coord_y)-1])
        elif d == Directions.R:
            coord_x.append(coord_x[len(coord_x)-1] + 1)
            coord_y.append(coord_y[len(coord_y)-1])
        ind_to_coord[i] = (coord_x[len(coord_x)-1], coord_y[len(coord_y) - 1])
        coord_to_ind[(coord_x[len(coord_x)-1], coord_y[len(coord_y) - 1])] = i
        i += 1
    return ind_to_coord, coord_to_ind, coord_x, coord_y 


def algorithmA(p: list[int]) -> list[Directions]:
    dir = []
    paux = prot(p, False)
    pb = prot_block(Block_type.SEP)
    pb.add_amino(0)
    for _ in range(0, paux.getBlocks()[0].getSize()):
        dir.append(Directions.D)
    P1, P2 = subroutine1(prot([pb] + paux.getBlocks()[1:-1] + [pb]))
    dir = dir + P1.fold(Block_type.Y_BLOCK)
    dir.append(Directions.R)
    dir = dir + P2.fold(Block_type.X_BLOCK, True)
    for _ in range(0, paux.getBlocks()[-1].getSize()):
        dir.append(Directions.U)
    return dir
    
str_seq = '0101000111100100000001101011000101101001'
seq = [1, 1, 1, 1, 3, 4, 2, 1, 7, 2, 1, 1, 1, 2, 3, 1, 1, 2, 1, 1, 2, 1, 0]
res = algorithmA(seq)
_, _, coord_x, coord_y = prot_coord2D(res)
color = []
for ch in str_seq:
    if ch == '1':
        color.append('black')
    else:
        color.append('white')
plt.plot(coord_x, coord_y, '-', c = 'black', zorder = 1)
plt.axis('equal')
plt.scatter(coord_x, coord_y, c = color, s = 100, edgecolors = 'black', zorder = 2)
plt.grid(color = 'black', linewidth=0.5)
plt.show()

str_seq = '010111000101000100101011'
seq = format_seq(str_seq)
print(seq)
res = algorithmA(seq)
_, _, coord_x, coord_y = prot_coord2D(res)
color = []
for ch in str_seq:
    if ch == '1':
        color.append('black')
    else:
        color.append('white')
plt.plot(coord_x, coord_y, '-', c = 'black', zorder = 1)
plt.axis('equal')
plt.scatter(coord_x, coord_y, c = color, s = 100, edgecolors = 'black', zorder = 2)
plt.grid(color = 'black', linewidth=0.5)
plt.show()