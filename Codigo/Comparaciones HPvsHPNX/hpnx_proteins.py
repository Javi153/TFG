from enum import Enum
import proteins as pts

#Directions = Enum('Directions', ['U', 'D', 'L', 'R', 'F', 'B'])
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
                
    def getType(self):
        return self._type
    
    def getN(self):
        return self._num_h
    
    def getSize(self):
        return self._num_aminos
    
    def setType(self, type1):
        self._type = type1

    def centrals(self):
        res = [0]
        acum = 1
        for i in range(1, len(self._seq)):
            if i % 2 == 0:
                res.append(acum)
            acum += self._seq[i]
        return res
    
class prot:
    _blocks : list[prot_block]
    _nx = 0
    _ny = 0
    _type = Block_type.Y_BLOCK
        
    def __init__(self, blocks, isBlock = True, type = Block_type.Y_BLOCK):
        self._nx = 0
        self._ny = 0
        self._type = type
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
        if self._blocks[-1].getN() == 1:
            self._blocks.pop()
        else:
            self._blocks[-1].del_amino()

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
            while i < len(seq) and ((i % 2 != 0 and seq[i] == 1) or (i % 2 == 0 and seq[i] % 2 != 0 and i != len(seq) - 1)):
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
            else:
                self._blocks.append(pb)
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
    
    def getType(self):
        return self._type
    
    def setType(self, type):
        self._type = type
    
    def centrals(self):
        acum = 0
        centrals = []
        for bl in self._blocks:
            if bl.getType() == self._type:
                centrals = centrals + list(map(lambda num: num + acum, bl.centrals()))
            acum += bl.getSize()
        return centrals