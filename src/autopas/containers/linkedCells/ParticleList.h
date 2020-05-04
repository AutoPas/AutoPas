//
// Created by lunaticcoding on 04.05.20.
//

#pragma AUTOPAS_PARTICLELIST_H

#include <vector>

template <class Type>
class ParticleList{

    using particleListImpType = std::vector<Type>;

    using iterator = typename particleListImpType::iterator;
    using const_iterator = typename particleListImpType::const_iterator;

    public:
        ParticleList<Type>() {
            dirty = false;
            particleListImp = std::vector<Type>();
        }

        Type* getReference(int index) {
            return &particleListImp[index];
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

        int size() {
            return particleListImp.size();
        }

        iterator begin() { return particleListImp.begin(); }
        iterator end() { return particleListImp.end(); }
        const_iterator begin() const { return particleListImp.begin(); }
        const_iterator end() const { return particleListImp.end(); }
        const_iterator cbegin() const { return particleListImp.cbegin(); }
        const_iterator cend() const { return particleListImp.cend(); }

    private:
        bool dirty;
        autopas::AutoPasLock particleListLock;
        std::vector<Type> particleListImp;
};

