#ifndef HP_PROT
#define HP_PROT
#include <iostream>
#include <vector>
#include <stack>

using namespace std;

enum block_type {X_BLOCK, Y_BLOCK, SEP};

enum direction {UP, DOWN, LEFT, RIGHT, FORWARD, BACKWARD};

class prot_block{
    private:
        block_type _type;
        vector<int> _seq;
        vector<direction> _dir;
        int _num1;
        int _num_amino;
    public:
        prot_block(){
            _type = SEP;
            _seq = vector<int>();
            _num1 = 0;
            _num_amino = 0;
            _dir = vector<direction>();
        }

        prot_block(const block_type& c_type, const vector<int>& c_seq) : _type(c_type), _seq(c_seq) {
            _num1 = 0;
            _num_amino = 0;
            for(int i = 0; i < c_seq.size(); ++i){
                if(i % 2 != 0){
                    _num1 += c_seq[i];
                }
                _num_amino += c_seq[i];
            }
            _dir = vector<direction>();
        }

        prot_block(const block_type& c_type) : _type(c_type){
            _seq = vector<int>();
            _num1 = 0;
            _num_amino = 0;
            _dir = vector<direction>();
        }

        prot_block(const prot_block& pb){
            _type = pb._type;
            _seq = vector<int>(pb._seq);
            _num1 = pb._num1;
            _num_amino = pb._num_amino;
            _dir = vector<direction>();
        }

        void add_seq(int am){
            if(_seq.size() % 2 != 0){
                _num1 += am;
            }
            _num_amino += am;
            _seq.push_back(am);
        }

        void del_amino(){
            if(_seq.size() != 0){
                if(_seq.size() % 2 == 0){
                    _num1 -= _seq[_seq.size() - 1];
                }
                _num_amino -= _seq[_seq.size() - 1];
                _seq.pop_back();
            }
        }

        block_type getType() const{
            return _type;
        }

        const vector<int>& getSeq() const{
            return _seq;
        }

        int getN() const{
            return _num1;
        }

        int size() const{
            return _num_amino;
        }

        const vector<direction>& fold_x(int size = 0, bool inverse = false){
            direction vert, hor1, hor2;
            if(inverse){
                vert = UP;
                hor1 = RIGHT;
                hor2 = LEFT;
            }
            else{
                vert = DOWN;
                hor1 = LEFT;
                hor2 = RIGHT;
            }
            if(_type == X_BLOCK){
                for(int i = 0; i < _seq.size(); ++i){
                    if(i % 2 == 0){
                        _dir.push_back(vert);
                    }
                    else{
                        for(int j = 0; j < _seq[i] / 2; ++j){
                            _dir.push_back(hor1);
                        }
                        _dir.push_back(vert);
                        for(int j = 0; j < _seq[i] / 2; ++j){
                            _dir.push_back(hor2);
                        }
                    }
                }
            }
        }

        const vector<direction>& fold_y(){

        }
};

class protein{
    private:
        vector<int> _seq;
        vector<prot_block> _blocks;
        bool _partitioned;
    public:
        protein(){
            _seq = vector<int>();
            _partitioned = false;
        }

        protein(const vector<int>& c_seq) : _seq(c_seq) {
            _partitioned = false;
        }

        protein(const protein& pr){
            _seq = vector<int>(pr._seq);
            _blocks = vector<prot_block>(pr._blocks);
            _partitioned = pr.isPartitioned();
        }

        protein(const vector<prot_block>& c_blocks) : _blocks(c_blocks){
            _seq = vector<int>();
            _partitioned = true;
        }

        void add_amino(int am){
            _seq.push_back(am);
            _partitioned = false;
        }

        void del_amino(){
            if(_seq.size() != 0){
                _seq.pop_back();
                _partitioned = false;
            }
        }

        void update_partition(){
            bool next_y = true;
            bool carry = false;
            block_type bt = Y_BLOCK;
            prot_block pb = prot_block(SEP);
            pb.add_seq(_seq[0]);
            _blocks = vector<prot_block>();
            _blocks.push_back(pb);
            int i = 1; //Suponemos que la primera ristra son ceros aunque este vacio
            while(i < _seq.size()){
                next_y ? bt = Y_BLOCK : bt = X_BLOCK;
                next_y = !next_y;
                pb = prot_block(bt);
                if(carry){
                    carry = !carry;
                    pb.add_seq(1);
                    ++i;
                }
                while(i < _seq.size() && ((i % 2 != 0 && _seq[i] == 1) || (i % 2 == 0 && _seq[i] % 2 != 0))){
                    pb.add_seq(_seq[i]);
                    ++i;
                }
                if(i % 2 == 0){
                    _blocks.push_back(pb);
                    if(i != _seq.size() - 1){
                        pb = prot_block(SEP);
                        pb.add_seq(_seq[i]);
                        _blocks.push_back(pb);
                        ++i;
                    }
                }
                else{
                    pb.add_seq(1);
                    _blocks.push_back(pb);
                    for(int j = 0; j < _seq[i] - 1; ++j){
                        next_y ? bt = Y_BLOCK : bt = X_BLOCK;
                        next_y = !next_y;
                        pb = prot_block(bt);
                        pb.add_seq(1);
                        _blocks.push_back(pb);
                    }
                    carry = true;
                }
            }
            _partitioned = true;
        }

        int Nx() const{
            if(_partitioned){
                int nx = 0;
                for(const prot_block& pb : _blocks){
                    if(pb.getType() == X_BLOCK){
                        nx += pb.getN();
                    }
                }
                return nx;
            }
            else{
                return 0;
            }
        }

        int Ny() const{
            if(_partitioned){
                int ny = 0;
                for(const prot_block& pb : _blocks){
                    if(pb.getType() == Y_BLOCK){
                        ny += pb.getN();
                    }
                }
                return ny;
            }
            else{
                return 0;
            }
        }

        int getSeqSize() const{
            return _seq.size();
        }

        int getBlockSize() const{
            return _blocks.size();
        }

        vector<prot_block>& getBlocks(int ini, int fin) const{
            vector<prot_block> vpb;
            int i = 0, num_blocks = 0;
            while(i < _blocks.size() && num_blocks < ini){
                if(_blocks[i].getType() != SEP){
                    num_blocks++;
                }
                ++i;
            }
            while(i < _blocks.size() && num_blocks < fin){
                if(_blocks[i].getType() != SEP){
                    num_blocks++;
                }
                vpb.push_back(_blocks[i]);
                ++i;
            }
            if(_blocks[i].getType() == SEP){
                vpb.push_back(_blocks[i]);
            }
            return vpb;
        }

        vector<prot_block>& getBlocks(int ini) const{
            vector<prot_block> vpb;
            int i = 0, num_blocks = 0;
            while(i < _blocks.size() && num_blocks < ini){
                if(_blocks[i].getType() != SEP){
                    num_blocks++;
                }
                ++i;
            }
            while(i < _blocks.size()){
                vpb.push_back(_blocks[i]);
                ++i;
            }
            return vpb;
        }

        int getMaxBlocks() const{
            return _blocks.size();
        }

        bool isPartitioned() const{
            return _partitioned;
        }
};
#endif