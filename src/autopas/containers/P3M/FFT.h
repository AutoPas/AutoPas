#pragma once

#include <complex>
#include <vector>
#include <math.h>
#include <array>

#include "autopas/utils/WrapOpenMP.h"
#include "autopas/utils/ThreeDimensionalMapping.h"

namespace autopas {
    
    class FFT {

        public:
        using RealGridType = typename std::vector<double>;
        using ComplexGridType = std::vector<std::complex<double>>;

        private:
        std::vector<std::complex<double>> buffer;
        std::vector<std::complex<double>> ws;
        ComplexGridType fftBuffer3D;
        bool correct_scaling;

        void rearange(std::vector<std::complex<double>> &x, unsigned int N){
            unsigned int bitReverse = N >> 1;
            for(unsigned int i = 1; i < N; i++){
                if(i < bitReverse){
                    swap(x[i], x[bitReverse]);
                }

                // add +1 reversed
                unsigned int mask = N >> 1;
                while(bitReverse & mask){
                    //sets leading 1 to 0
                    bitReverse &= ~mask;
                    mask = mask >> 1;
                }
                bitReverse |= mask;
            }
        }

        // requires x = result in initial call
        void fft(std::vector<std::complex<double>> &x, unsigned int N, bool forward){
            rearange(x, N);

            /*for(unsigned int i = 0; i < N; i++){
                std::cout << x[i] << ", ";
            }
            std::cout << std::endl;*/

            unsigned int period = 2;

            std::complex<double> w_n;
            

            while(period <= N){
                if(forward){
                    w_n = std::exp(std::complex(0.,-2*M_PI/ period));
                }else{
                    w_n = std::exp(std::complex(0.,2*M_PI/ period));
                }
                

                for(unsigned int i = 0; i < N/period; i++){
                    std::complex<double> w = 1;
                    for(unsigned int k = 0; k < period/2; k++){
                        auto tmp = w * x[i*period + k + (period/2)];
                        x[i*period + k + (period/2)] = x[i*period + k] - tmp;
                        x[i*period + k] += tmp;

                        w *= w_n;
                    }
                }
                period *= 2;
            }
        }

        public:

        FFT(){}

        FFT(std::array<unsigned int, 3> &N, bool correct_scaling = false) : correct_scaling(correct_scaling){
            unsigned int maxGridDim = 0;
            for (int i = 0; i < 3; i++){
                //power of two check
                if(N[i] & (N[i]-1) || N[i] < 1){
                    //TODO throw error
                    return;
                }
                if(N[i] > maxGridDim){
                    maxGridDim = N[i];
                }
            }
            buffer = std::vector<std::complex<double>>(maxGridDim);

            fftBuffer3D = ComplexGridType(N[0] * N[1] * N[2]);
            
        }

        //TOOD move construction

        void forward(std::vector<std::complex<double>> &x, unsigned int N){
            fft(x, N, true);
            if(correct_scaling){
                for(unsigned int i = 0; i < N; i++){
                    x[i] /= std::sqrt(2*M_PI);
                }
            }
        }

        void backward(std::vector<std::complex<double>> &x, unsigned int N){
            fft(x, N, false);
            if(correct_scaling){
                for(unsigned int i = 0; i < N; i++){
                    x[i] = x[i] * std::sqrt(2*M_PI) / std::complex<double>(N);
                }
            }else{
                for(unsigned int i = 0; i < N; i++){
                    x[i] = x[i] / std::complex<double>(N);
                }
            }
            
        }

        /*void forward3DdirectionAt(ComplexGridType &in, ComplexGridType &out, std::array<unsigned int, 3> &grid_dims, int direction, int xIndx, int yIndx, int zIndx){
            for(unsigned int m = 0; m < grid_dims[direction]; m++){
                // copy values into buffer
                switch(direction){
                    case 0:
                        buffer[m] = in[m][yIndx][zIndx];
                        break;
                    case 1:
                        buffer[m] = in[xIndx][m][zIndx];
                        break;
                    case 2:
                        buffer[m] = in[xIndx][yIndx][m];
                        break;
                }
                        
            }
            forward(buffer, grid_dims[direction], outBuffer);
            for(unsigned int m = 0; m < grid_dims[direction]; m++){
                // copy values into buffer
                switch(direction){
                    case 0:
                        out[m][yIndx][zIndx] = outBuffer[m];
                        break;
                    case 1:
                        out[xIndx][m][zIndx] = outBuffer[m];
                        break;
                    case 2:
                        out[xIndx][zIndx][m] = outBuffer[m];
                        break;
                }
            }
        }*/

        void forward3D(RealGridType &in, ComplexGridType &out, std::array<unsigned int, 3> &grid_dims){
            AUTOPAS_OPENMP(parallel for schedule(static))
            for(unsigned int i = 0; i < grid_dims[0] * grid_dims[1] * grid_dims[2]; i++){
                fftBuffer3D[i] = std::complex<double>(in[i]);
            }


            AUTOPAS_OPENMP(parallel for schedule(static) firstprivate(buffer))
            for (unsigned int z = 0; z < grid_dims[2]; z++){
                // do all FFTs for one x-Plane
                for(unsigned int y = 0; y < grid_dims[1]; y++){
                    // FFTs in x-direction
                    for(unsigned int x = 0; x < grid_dims[0]; x++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        buffer[x] = fftBuffer3D[index1d];
                    }
                    forward(buffer, grid_dims[0]);
                    for(unsigned int x = 0; x < grid_dims[0]; x++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        fftBuffer3D[index1d] = buffer[x];
                    }
                }
                for(unsigned int x = 0; x < grid_dims[0]; x++){
                    // FFTs in y-direction on the already transformed values
                    for(unsigned int y = 0; y < grid_dims[1]; y++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        buffer[y] = fftBuffer3D[index1d];
                    }
                    forward(buffer, grid_dims[1]);
                    for(unsigned int y = 0; y < grid_dims[1]; y++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        fftBuffer3D[index1d] = buffer[y];
                    }
                }
            }

            //FFTs in z-direction
            AUTOPAS_OPENMP(parallel for schedule(static) firstprivate(buffer))
            for(unsigned y = 0; y < grid_dims[1]; y++){
                for(unsigned x = 0; x < grid_dims[0]; x++){
                    for(unsigned z = 0; z < grid_dims[2]; z++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        buffer[z] = fftBuffer3D[index1d];
                    }
                    forward(buffer, grid_dims[2]);
                    for(unsigned z = 0; z < grid_dims[2]; z++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        out[index1d] = buffer[z];
                    }
                }
            }
        }

        /*void backward3DdirectionAt(ComplexGridType &in, ComplexGridType &out, std::array<unsigned int, 3> &grid_dims, int direction, int xIndx, int yIndx, int zIndx){
            for(unsigned int m = 0; m < grid_dims[direction]; m++){
                // copy values into buffer
                switch(direction){
                    case 0:
                        buffer[m] = in[m][yIndx][zIndx];
                        break;
                    case 1:
                        buffer[m] = in[xIndx][m][zIndx];
                        break;
                    case 2:
                        buffer[m] = in[xIndx][yIndx][m];
                        break;
                }
            }
            backward(buffer, grid_dims[direction], outBuffer);
            for(unsigned int m = 0; m < grid_dims[direction]; m++){
                // copy values into buffer
                switch(direction){
                    case 0:
                        out[m][yIndx][zIndx] = outBuffer[m];
                        break;
                    case 1:
                        out[xIndx][m][zIndx] = outBuffer[m];
                        break;
                    case 2:
                        out[xIndx][yIndx][m] = outBuffer[m];
                        break;
                }
            }
        }*/

        void backward3D(ComplexGridType &in, RealGridType &out, std::array<unsigned int, 3> &grid_dims){
            AUTOPAS_OPENMP(parallel for schedule(static) firstprivate(buffer))
            for (unsigned int z = 0; z < grid_dims[2]; z++){
                // do all FFTs for one x-Plane
                for(unsigned int y = 0; y < grid_dims[1]; y++){
                    // FFTs in x-direction
                    for(unsigned int x = 0; x < grid_dims[0]; x++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        buffer[x] = in[index1d];
                    }
                    backward(buffer, grid_dims[0]);
                    for(unsigned int x = 0; x < grid_dims[0]; x++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        fftBuffer3D[index1d] = buffer[x];
                    }
                }
                for(unsigned int x = 0; x < grid_dims[0]; x++){
                    // FFTs in y-direction on the already transformed values
                    for(unsigned int y = 0; y < grid_dims[1]; y++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        buffer[y] = fftBuffer3D[index1d];
                    }
                    backward(buffer, grid_dims[1]);
                    for(unsigned int y = 0; y < grid_dims[1]; y++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        fftBuffer3D[index1d] = buffer[y];
                    }
                }
            }

            //FFTs in z-direction
            AUTOPAS_OPENMP(parallel for schedule(static) firstprivate(buffer))
            for(unsigned y = 0; y < grid_dims[1]; y++){
                for(unsigned x = 0; x < grid_dims[0]; x++){
                    for(unsigned z = 0; z < grid_dims[2]; z++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        buffer[z] = fftBuffer3D[index1d];
                    }
                    backward(buffer, grid_dims[2]);
                    for(unsigned z = 0; z < grid_dims[2]; z++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        fftBuffer3D[index1d] = buffer[z];
                    }
                }
            }

            // Make output real
            AUTOPAS_OPENMP(parallel for schedule(static))
            for(unsigned int i = 0; i < grid_dims[0] * grid_dims[1] * grid_dims[2]; i++){
                out[i] = std::real(fftBuffer3D[i]);
            }
        }

        
    };
}