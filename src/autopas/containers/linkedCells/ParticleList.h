//
// Created by lunaticcoding on 04.05.20.
//

#pragma AUTOPAS_PARTICLELIST_H

#include <vector>

template <class Type>
class ParticleList{
    public:
        ParticleList<Type>() {
            dirty = false;
        }

        ParticleList<Type>(int n) {
            dirty = false;
            particleListImp = std::vector<Type>(n);
        }

        Type get(int index) {
            return particleListImp[index];
        }

        void set(int index, Type value) {
            particleListImp[index] = value;
        }

        void push_back(Type &value) {
            dirty &= particleListImp.capacity() == particleListImp.size();
            particleListLock.lock();
            particleListImp.push_back(value);
            particleListLock.unlock();
        }

        Type pop_back() {
            return particleListImp.pop_back();
        }

        Type* begin(){
            return particleListImp.begin();
        }

        Type* end(){
            return particleListImp.end();
        }

    private:
        bool dirty;
        autopas::AutoPasLock particleListLock;
        std::vector<Type> particleListImp;
};

