#pragma once
#include <iostream>
#include <vector>

struct Hydrogen {
    int id;
    std::vector<double> x;
    std::vector<double> v;

    // ParticleBase() : _r({0.0, 0.0, 0.0}), _v({0., 0., 0.}), _f({0.0, 0.0, 0.0}), _id(0) {}

    Hydrogen() : id(0), x({1.0, 2.0, 3.0}), v({1.0,3.0,5.0}){
        std::cout << "hydrogen constructor with no arguments\n";
    }
    Hydrogen(int id_, std::vector<double> x_, std::vector<double> v_) : id(id_), x({x_}), v({v_}) {
        std::cout << "hydrogen constructor with arguments\n";
    }

    void update_x() {
        for(auto i = 0; i < x.size(); i++) {
            x.at(i) = x.at(i) + 5;
        }
    }
    void update_v() {
        for(auto i = 0; i < v.size(); i++) {
            v.at(i) = v.at(i) * 2;
        }
    }
};