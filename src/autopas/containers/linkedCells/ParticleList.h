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

        bool isDirty() {
            return dirty;
        }

        void markAsClean() {
            dirty = false;
        }

        void set(int index, Type value) {
            particleListImp[index] = value;
        }

        void push_back(Type &value) {
            particleListLock.lock();
            dirty &= particleListImp.capacity() == particleListImp.size();
            particleListImp.push_back(value);
            particleListLock.unlock();
        }

        // TODO push back all or pass vector/template  push_back({p1, p2, p3}); - later
//        void push_back(Type ...&value) {
//            particleListLock.lock();
//            dirty &= particleListImp.capacity() == particleListImp.size();
//            particleListImp.push_back(value);
//            particleListLock.unlock();
//        }

        void emblace_back(Type &value) {
            // TODO factor out into function
            particleListLock.lock();
            dirty &= particleListImp.capacity() == particleListImp.size();
            particleListImp.push_back(value);
            particleListLock.unlock();
        }

        Type pop_back(Type &value) {
            particleListLock.lock();
            dirty &= particleListImp.capacity() == particleListImp.size();
            particleListImp.pop_back(value);
            particleListLock.unlock();
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
    // TODO index since markAsClean() / maybe not now
        bool dirty;
        autopas::AutoPasLock particleListLock;
        std::vector<Type> particleListImp;
};

