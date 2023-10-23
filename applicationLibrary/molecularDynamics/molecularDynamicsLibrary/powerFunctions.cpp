//
// Created by kay on 19.10.23.
//

// a class for testing
#include <iostream>
#include <utility>
#include <immintrin.h>
#include <vector>
#include <map>

#ifndef AUTOPAS_POWERFUNCTIONS_H
#define AUTOPAS_POWERFUNCTIONS_H
struct AdditionChain{
  std::vector<unsigned> additionChain;
};

// best case is that the result b^m and b^n is already computed in the loop

double additionExpo(const AdditionChain& chain, double b){
  std::map<unsigned, double> computedChain{};
  computedChain[1]=b;
  for( int j = 1; j < chain.additionChain.size();j++){
    auto k = chain.additionChain[j] - chain.additionChain[j-1];
    std::cout << "k is:" << k << std::endl;
    computedChain[chain.additionChain[j]]=computedChain[chain.additionChain[j-1]]* (computedChain[k]);
    std::cout << "midvalue for j = " << j << " is: " << computedChain[chain.additionChain[j]] << std::endl;
  }
  std::cout << "endvalue is: " << computedChain[chain.additionChain.back()] << std::endl;
  return computedChain[computedChain.size()];
  }

// exponentationBySquaring
void exponentationBySquaringSIMD(unsigned m, unsigned n, const double * d) {
  __m256d squared;
  __m256d result_m {_mm256_set1_pd(1.)};
  __m256d result_n{_mm256_set1_pd(1.)};
  __m256d b = _mm256_load_pd(d);
    if(m & 1 || n & 1){
      squared = _mm256_sqrt_pd(b);
    }
    //start from squared
    m<<=1;
    n<<=1;

  while (m > 0) {
    if (m & 1) {
      result_m = _mm256_mul_pd(result_m,b) ;
    }
    b = _mm256_mul_pd(b,b) ;
    m >>= 1;
  }

   b = _mm256_load_pd(d);
  while (n > 0) {
    if (n & 1) {
      result_n = _mm256_mul_pd(result_n,b) ;
    }
    b = _mm256_mul_pd(b,b) ;
    n >>= 1;
  }
  double d1[4];
  double d2[4];
  _mm256_storeu_pd(d1,result_n);
  _mm256_storeu_pd(d2,result_m);
  for (int i = 0; i < 4; i++) {
    std::cout << "result[" << i << "] = " << "result_n: "<< d1[i] <<"result_m: " << d2[i] <<std::endl;
  }


}

// exponentationBySquaring
void exponentationBySquaringSIMD2(unsigned m, unsigned n, const double * d) {
  __m256d squared;
  __m256d result_m {_mm256_set1_pd(1.)};
  __m256d result_n{_mm256_set1_pd(1.)};
  __m256d b = _mm256_load_pd(d);
  if(m & 1 || n & 1){
    squared = _mm256_sqrt_pd(b);
  }
  //start from squared
  m<<=1;
  n<<=1;

  while (m > 0 || n > 0) {
    if (m & 1) {
      result_m = _mm256_mul_pd(result_m,b) ;
    }
    if (m & 1){
      result_n = _mm256_mul_pd(result_m,b) ;

    }
    b = _mm256_mul_pd(b,b) ;
    m >>= 1;
    n >>=1;
  }


  double d1[4];
  double d2[4];
  _mm256_storeu_pd(d1,result_n);
  _mm256_storeu_pd(d2,result_m);
  for (int i = 0; i < 4; i++) {
    std::cout << "result[" << i << "] = " << "result_n: "<< d1[i] <<"result_m: " << d2[i] <<std::endl;
  }


}



std::pair<double ,double > exponentationBySquaring(size_t m, size_t n, double base){

  double result_m = 1.0;
  double result_n = 1.0;
  double b = base;

  while (m > 0) {
    if (m & 1) {
      result_m *= b;
    }
    b *= b;
    m >>= 1;
  }

  b=base;
  while (n > 0) {
    if (n & 1) {
      result_n *= b;
    }
    b *= b;
    n >>= 1;
  }
  return std::make_pair(result_m,result_n);

}
// probably a bit faster, uses one loop
std::pair<double ,double > exponentationBySquaring2(size_t m, size_t n, double base){

  double result_m = 1.0;
  double result_n = 1.0;
  double b = base;
  // option 1 compute bigger and save smaller when computing
  // option 2 compute smaller and then compute bigger

  while (m > 0) {
    if (m & 1) {
      result_m *= b;
    }
    b *= b;
    m >>= 1;
  }

  b=base;
  while (n > 0) {
    if (n & 1) {
      result_n *= b;
    }
    b *= b;
    n >>= 1;
  }
  return std::make_pair(result_m,result_n);

}

int main(){
  /*
  std::vector<unsigned> chain{1,2,3,6,9,12,15};
  AdditionChain additionChain{chain};
  additionExpo(additionChain,3.125);
  return 0;
   */
 const double base [4] = {1,2,3,2.125};
  exponentationBySquaringSIMD(12,5, base);
  exponentationBySquaringSIMD2(12,5, base);

  return 1;
  auto [a,b] = exponentationBySquaring(12,5,2.4);
    std::cout << a << "  " << b << std::endl;

}

#endif  // AUTOPAS_POWERFUNCTIONS_H
