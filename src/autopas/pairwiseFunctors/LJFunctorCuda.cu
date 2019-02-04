/**
 * @file LJFunctorCuda.h
 *
 * @date 9.1.2019
 * @author jspahl
 */
#include "autopas/pairwiseFunctors/LJFunctorCuda.cuh"

__constant__ constants global_constants;

__device__
double3 bodyBodyF(double3 i,
		double3 j, double3 fi) {
  double drx = i.x - j.x;
  double dry = i.y - j.y;
  double drz = i.z - j.z;

  double dr2 = drx* drx + dry * dry + drz* drz;

  if (dr2 > global_constants.cutoffsquare | dr2 == 0.0) {
  	return fi;
  }

  double invdr2 = 1. / dr2;
  double lj6 = global_constants.sigmasquare * invdr2;
  lj6 = lj6 * lj6 * lj6;
  double lj12 = lj6 * lj6;
  double lj12m6 = lj12 - lj6;
  double fac = global_constants.epsilon24 * (lj12 + lj12m6) * invdr2;

  fi.x += drx * fac;
  fi.y += dry * fac;
  fi.z += drz * fac;

  return fi;
}


template <int block_size>
__global__
void SoAFunctorNoN3(int N, double* posX, double* posY, double* posZ, double* forceX, double* forceY, double* forceZ){
	 __shared__ double3 block_pos[block_size];
	 int i, tile;
	 int tid = blockIdx.x * blockDim.x + threadIdx.x;
	 if(tid >= N){
		 return;
	 }

	 double3 myposition;
	 myposition.x = posX[tid];
	 myposition.y = posY[tid];
	 myposition.z = posZ[tid];
	 double3 myf = {forceX[tid], forceY[tid], forceZ[tid]};

	 for(i = 0, tile = 0; i < N; i+=block_size, ++tile){
		 int idx = tile * blockDim.x + threadIdx.x;

		 block_pos[threadIdx.x] = {posX[idx], posY[idx], posZ[idx]};
		 __syncthreads();

		  for(int j = 0; j < blockDim.x; ++j){
			  myf = bodyBodyF(myposition, block_pos[j], myf);
		  }
		 __syncthreads();
	 }

	 forceX[tid] = myf.x;
	 forceY[tid] = myf.y;
	 forceZ[tid] = myf.z;
}

template <int block_size>
__global__
void SoAFunctorNoN3Pair(int N, double* posX, double* posY, double* posZ, double* forceX, double* forceY, double* forceZ,
		int M, double* posX2, double* posY2, double* posZ2){
	 __shared__ double3 block_pos[block_size];
	 int i, tile;
	 int tid = blockIdx.x * blockDim.x + threadIdx.x;
	 if(tid >= N){
		 return;
	 }

	 double3 myposition;
	 myposition.x = posX[tid];
	 myposition.y = posY[tid];
	 myposition.z = posZ[tid];
	 double3 myf = {forceX[tid], forceY[tid], forceZ[tid]};

	 for(i = 0, tile = 0; i < M; i+=block_size, ++tile){
		 int idx = tile * blockDim.x + threadIdx.x;

		 block_pos[threadIdx.x] = {posX2[idx], posY2[idx], posZ2[idx]};
		 __syncthreads();

		  for(int j = 0; j < blockDim.x; ++j){
			  myf = bodyBodyF(myposition, block_pos[j], myf);
		  }
		 __syncthreads();
	 }

	 forceX[tid] = myf.x;
	 forceY[tid] = myf.y;
	 forceZ[tid] = myf.z;
}

template <int block_size>
__global__
void AoSFunctorNoN3(int N, double* particles){
	 __shared__ double3 block_pos[block_size];
	 int i, tile;
	 int tid = blockIdx.x * blockDim.x + threadIdx.x;
	 if(tid >= N){
		 return;
	 }

	 double3 myposition;
	 myposition.x = particles[6 * tid + 0];
	 myposition.y = particles[6 * tid + 1];
	 myposition.z = particles[6 * tid + 2];
	 double3 myf = {particles[6 * tid + 3],particles[6 * tid + 4],particles[6 * tid + 5]};

	 for(i = 0, tile = 0; i < N; i+=block_size, ++tile){
		 int idx = tile * blockDim.x + threadIdx.x;

		 block_pos[threadIdx.x] = {particles[6 * idx + 0], particles[6 * idx + 1], particles[6 * idx + 2]};
		 __syncthreads();

		  for(int j = 0; j < blockDim.x; ++j){
			  myf = bodyBodyF(myposition, block_pos[j], myf);
		  }
		 __syncthreads();
	 }

	 particles[6 * tid + 3] = myf.x;
	 particles[6 * tid + 4] = myf.y;
	 particles[6 * tid + 5] = myf.z;
}

template <int block_size>
__global__
void AoSFunctorNoN3Pair(int N, int M, double* particles1, double* particles2){
	 int i, tile;
	 int tid = blockIdx.x * blockDim.x + threadIdx.x;
	 __shared__ double3 block_pos2[32];

	 double3 myposition;
	 myposition.x = particles1[6 * tid + 0];
	 myposition.y = particles1[6 * tid + 1];
	 myposition.z = particles1[6 * tid + 2];
	 double3 myf = {0,0,0 };

	 for(i = 0, tile = 0; i < M; i+=block_size, ++tile){
		 int idx = tile * blockDim.x + threadIdx.x;

		 block_pos2[threadIdx.x] = {particles2[6 * idx + 0], particles2[6 * idx + 1], particles2[6 * idx + 2]};
		 __syncthreads();

		  for(int j = 0; j < blockDim.x; ++j){
			  myf = bodyBodyF(myposition, block_pos2[j], myf);
		  }
		 __syncthreads();
	 }

	 particles1[6 * tid + 3] = myf.x;
	 particles1[6 * tid + 4] = myf.y;
	 particles1[6 * tid + 5] = myf.z;
}

void AoSFunctorNoN3Wrapper(int N, double* particles){
	AoSFunctorNoN3<32><<<N/32 + 1,32>>>(N, particles);
}

void AoSFunctorNoN3PairWrapper(int N, int M, double* particles1, double* particles2){
	AoSFunctorNoN3Pair<32><<<N/32 + 1,32>>>(N, M, particles1, particles2);
}

void SoAFunctorNoN3Wrapper(int N, double* posX, double* posY, double* posZ, double* forceX, double* forceY, double* forceZ){
	SoAFunctorNoN3<32><<<N/32 + 1,32>>>(N, posX, posY, posZ, forceX, forceY, forceZ);
}

void SoAFunctorNoN3PairWrapper(int N, double* posX, double* posY, double* posZ, double* forceX, double* forceY, double* forceZ,
		int M, double* posX2, double* posY2, double* posZ2){
	SoAFunctorNoN3Pair<32><<<N/32 + 1,32>>>(N, posX, posY, posZ, forceX, forceY, forceZ, M, posX2, posY2, posZ2);
}

void loadConstants(double cutoffsquare, double epsilon24, double sigmasquare){

	constants c;
	c.cutoffsquare = cutoffsquare;
	c.epsilon24 = epsilon24;
	c.sigmasquare = sigmasquare;

	cudaMemcpyToSymbol(global_constants, &c, sizeof(constants));
}

//template __global__ void SoAFunctorNoN3<32>(int N, double* posX, double* posY, double* posZ, double* forceX, double* forceY, double* forceZ);

