#include "autopas/pairwiseFunctors/LJFunctorCuda.cuh"

__constant__ constants global_constants;

__device__
double3 bodyBodyF(double3 i,
		double3 j, double3 fi, bool newton3) {
  double drx = i.x - j.x;
  double dry = i.y - j.y;
  double drz = i.z - j.z;

  double dr2 = drx* drx + dry * dry + drz* drz;

  if (dr2 > global_constants.cutoffsquare) {
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

  if (newton3) {
    //j.subF(f);
  }
  return fi;
}

__device__
double3 tile_calculation(double3 myposition, double3 f, bool newton3){
	  int i;
	  __shared__ double3 all_pos[32];
	  for(i = 0; i < blockDim.x; ++i){
		  f = bodyBodyF(myposition, all_pos[i], f, newton3);
	  }
	  return f;
}

__global__
void AoSFunctorCuda(int N, int p, double3* rs, double3* fs, bool newton3){
	 extern __shared__ double3 all_pos[];
	 int i, tile;
	 int tid = blockIdx.x * blockDim.x + threadIdx.x;
	 double3 myposition = rs[tid];
	 double3 myf = fs[tid];
	 for(i = 0, tile = 0; i < N; i+=p, ++tile){
		 int idx = tile * blockDim.x + threadIdx.x;
		 all_pos[threadIdx.x] = rs[idx];
		 __syncthreads();
		 myf = tile_calculation(myposition, myf, newton3);
		 __syncthreads();
	 }

	 fs[tid] = myf;
}

__global__
void SoAFunctorNoN3(int N, double* posX, double* posY, double* posZ, double* forceX, double* forceY, double* forceZ){

}

__global__
void AoSFunctorNoN3(int N, double* particles){
	 extern __shared__ double3 all_pos[];
	 int i, tile;
	 int tid = blockIdx.x * blockDim.x + threadIdx.x;
	 double3 myposition;
	 myposition.x = particles[6 * tid + 0];
	 myposition.y = particles[6 * tid + 1];
	 myposition.z = particles[6 * tid + 2];

	 double3 myf = {particles[6 * tid + 3],0,0 };
	 for(i = 0, tile = 0; i < N; i+=32, ++tile){
		 int idx = tile * blockDim.x + threadIdx.x;
		 all_pos[threadIdx.x] = {particles[6 * idx + 0], particles[6 * idx + 1], particles[6 * idx + 2]};
		 __syncthreads();
		 myf = tile_calculation(myposition, myf, false);
		 __syncthreads();
	 }

	 particles[6 * tid + 3] = myf.x;
	 particles[6 * tid + 4] = myf.y;
	 particles[6 * tid + 5] = myf.z;
}

void AoSFunctorNoN3Wrapper(int N, double* particles){
	AoSFunctorNoN3<<<N/32 + 1,32, N>>>(N, particles);
}

void AoSFunctorNoN3Wrapper(int N, int M, double* particles1, double* particles2){
	//AoSFunctorNoN3<<<N/32 + 1,32, N>>>(N, particles1, particles2);
}

void SoAFunctorNoN3Wrapper(int N, double* posX, double* posY, double* posZ, double* forceX, double* forceY, double* forceZ){
	SoAFunctorNoN3<<<N/32 + 1,32>>>(N, posX, posY, posZ, forceX, forceY, forceZ);
}

void loadConstants(double cutoffsquare, double epsilon24, double sigmasquare){

	constants c;
	c.cutoffsquare = cutoffsquare;
	c.epsilon24 = epsilon24;
	c.sigmasquare = sigmasquare;

	cudaMemcpyToSymbol(global_constants, &c, sizeof(constants));
}


