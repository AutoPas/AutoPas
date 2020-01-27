#include <stdio.h>

int main(int argc, char **argv){
    cudaDeviceProp dP;

    int rc = cudaGetDeviceProperties(&dP, 0);
    if(rc != cudaSuccess) {
        cudaError_t error = cudaGetLastError();
        printf("CUDA error: %s", cudaGetErrorString(error));
        return rc; /* Failure */
    }

    printf("%d%d", dP.major, dP.minor);
    return 0;
}
 
