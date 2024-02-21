#include <iostream>
#include <vector>
#include <utility>
#include "HP_prot.h"

using namespace std;

bool cond(int a, int b, int A, int B){
    return a < b || (a == b && A < B);
}

int Mxy(const vector<prot_block>& p1, const vector<prot_block>& p2, int min = true){
    int p1x = 0, p2y = 0;
    for(const prot_block& pb : p1){
        if(pb.getType() == X_BLOCK){
            p1x += pb.getN();
        }
    }
    for(const prot_block& pb : p2){
        if(pb.getType() == Y_BLOCK){
            p2y += pb.getN();
        }
    }
    if(min){
        return std::min(p1x, p2y);
    }
    else{
        return std::max(p1x, p2y);
    }
}

int Myx(const vector<prot_block>& p1, const vector<prot_block>& p2, int min = true){
    return Mxy(p2, p1, min);
}

pair<protein, protein> subrutine1(const protein& pr){
    vector<prot_block> B1, B2, P1, P2;
    int e, E, m1, m2, M1, M2;
    if(pr.isPartitioned()){
        B1 = pr.getBlocks(0, 1);
        B2 = pr.getBlocks(1);
        m1 = Mxy(B1, B2);
        m2 = Myx(B1, B2);
        M1 = Mxy(B1, B2, false);
        M2 = Myx(B1, B2, false);
        if(m1 > m2){
            P1 = vector<prot_block>(B2);
            P2 = vector<prot_block>(B1);
            e = m1;
            E = M1;
        }
        else{
            P1 = vector<prot_block>(B1);
            P2 = vector<prot_block>(B2);
            e = m2;
            E = M2;
        }
        for(int i = 3; i < pr.getMaxBlocks(); ++i){
            B1 = pr.getBlocks(0, i-1);
            B2 = pr.getBlocks(i);
            m1 = Mxy(B1, B2);
            m2 = Myx(B1, B2);
            M1 = Mxy(B1, B2, false);
            M2 = Myx(B1, B2, false);
            if(cond(e, m1, E, M1)){
                P1 = vector<prot_block>(B2);
                P2 = vector<prot_block>(B1);
                e = m1;
                E = M1;
            }
            else{
                P1 = vector<prot_block>(B1);
                P2 = vector<prot_block>(B2);
                e = m2;
                E = M2;
            }
        }
    }
    return make_pair(protein(P1), protein(P2));
}

protein& algorithmA(protein& pr){
    pr.update_partition();
    pair<protein, protein> aux = subrutine1(pr);
    protein B1 = aux.first, B2 = aux.second;
    B1.fold_x()
}

void algorithmB(){}

void algorithmC(){}

int main(){
    return 0;
}