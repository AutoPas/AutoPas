/**
 * @file LJFunctorCuda.h
 *
 * @date 14.12.2018
 * @author jspahl
 */

#pragma once

#include "cuda_runtime.h"


struct constants{
  double cutoffsquare;
  double epsilon24;
  double sigmasquare;
};

void loadConstants(double cutoffsquare, double epsilon24, double sigmasquare);

void AoSFunctorNoN3Wrapper(int N, double* particles);
void AoSFunctorNoN3Wrapper(int N, int M, double* particles1, double* particles2);

void SoAFunctorNoN3Wrapper(int N, double* posX, double* posY, double* posZ, double* forceX, double* forceY, double* forceZ);
