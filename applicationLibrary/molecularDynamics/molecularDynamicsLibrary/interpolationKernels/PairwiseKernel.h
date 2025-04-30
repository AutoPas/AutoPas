/**
 * @file PairwiseKernel.h
 * 
 * @date 29.04.2025
 * @author Luis Gall
 */

#pragma once

namespace mdLib {

template <class CRTP_T>
class PairwiseKernel {

public:
    explicit PairwiseKernel() {};

    virtual double calculate(double dr2) = 0;
};


} // mdLib