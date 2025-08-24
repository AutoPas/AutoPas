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
        std::vector<std::complex<double>> inBuffer;
        std::vector<std::complex<double>> outBuffer;
        std::vector<std::complex<double>> ws;
        ComplexGridType fftBuffer3D;
        bool correct_scaling;

        // requires x = result in initial call
        void fft(std::vector<std::complex<double>> &x, unsigned int N, int start, std::vector<std::complex<double>> &result, bool forward, int stride = 1){
            if(N == 1){
                //result[start] = x[start];
                return;
            }
            //split
            /*for(unsigned int i = 0; i < N/2; i++){
                result[start + i] = x[start + 2*i];
                result[start + N/2 + i] = x[start + 2*i+1];
            }*/
            fft(result, N/2, start, x, forward, 2*stride);
            fft(result, N/2, start + stride, x, forward, 2*stride);

            std::complex<double> w_n;
            if(forward){
                w_n = std::exp(std::complex(0.,-2*M_PI/ N));
            }else{
                w_n = std::exp(std::complex(0.,2*M_PI/ N));
            }
            std::complex<double> w = std::complex(1.,0.);

            for(unsigned int k = 0; k < (N/2); k++){
                result[start + k*stride] = x[start + 2*k*stride] + w * x[start + (2*k + 1)*stride];
                result[start + (N/2 + k)*stride] = x[start + 2*k*stride] - w * x[start + (2*k + 1)*stride];
                w *=  w_n;
            }
        }

        public:

        FFT(){}

        FFT(std::array<unsigned int, 3> &N, bool correct_scaling = false) : correct_scaling(correct_scaling){
            unsigned int maxGridDim = 0;
            for (int i = 0; i < 3; i++){
                //powerOfTwo(N[i]);
                if(N[i] > maxGridDim){
                    maxGridDim = N[i];
                }
            }
            inBuffer = std::vector<std::complex<double>>(maxGridDim);
            outBuffer = std::vector<std::complex<double>>(maxGridDim);

            fftBuffer3D = ComplexGridType(N[0] * N[1] * N[2]);

            ws = std::vector<std::complex<double>>((int) std::log2(maxGridDim));
            
        }

        //TOOD move construction

        void forward(std::vector<std::complex<double>> &x, unsigned int N, std::vector<std::complex<double>> &result){
            result = x;
            fft(x, N, 0, result, true);
            if(correct_scaling){
                for(unsigned int i = 0; i < N; i++){
                    result[i] /= std::sqrt(2*M_PI);
                }
            }
        }

        void backward(std::vector<std::complex<double>> &x, unsigned int N, std::vector<std::complex<double>> &result){
            result = x;
            fft(x, N, 0, result, false);
            if(correct_scaling){
                for(unsigned int i = 0; i < N; i++){
                    result[i] = result[i] * std::sqrt(2*M_PI) / std::complex<double>(N);
                }
            }else{
                for(unsigned int i = 0; i < N; i++){
                    result[i] = result[i] / std::complex<double>(N);
                }
            }
            
        }

        /*void forward3DdirectionAt(ComplexGridType &in, ComplexGridType &out, std::array<unsigned int, 3> &grid_dims, int direction, int xIndx, int yIndx, int zIndx){
            for(unsigned int m = 0; m < grid_dims[direction]; m++){
                // copy values into buffer
                switch(direction){
                    case 0:
                        inBuffer[m] = in[m][yIndx][zIndx];
                        break;
                    case 1:
                        inBuffer[m] = in[xIndx][m][zIndx];
                        break;
                    case 2:
                        inBuffer[m] = in[xIndx][yIndx][m];
                        break;
                }
                        
            }
            forward(inBuffer, grid_dims[direction], outBuffer);
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


            AUTOPAS_OPENMP(parallel for schedule(static) firstprivate(inBuffer, outBuffer))
            for (unsigned int z = 0; z < grid_dims[2]; z++){
                // do all FFTs for one x-Plane
                for(unsigned int y = 0; y < grid_dims[1]; y++){
                    // FFTs in x-direction
                    for(unsigned int x = 0; x < grid_dims[0]; x++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        inBuffer[x] = fftBuffer3D[index1d];
                    }
                    forward(inBuffer, grid_dims[0], outBuffer);
                    for(unsigned int x = 0; x < grid_dims[0]; x++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        fftBuffer3D[index1d] = outBuffer[x];
                    }
                }
                for(unsigned int x = 0; x < grid_dims[0]; x++){
                    // FFTs in y-direction on the already transformed values
                    for(unsigned int y = 0; y < grid_dims[1]; y++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        inBuffer[y] = fftBuffer3D[index1d];
                    }
                    forward(inBuffer, grid_dims[1], outBuffer);
                    for(unsigned int y = 0; y < grid_dims[1]; y++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        fftBuffer3D[index1d] = outBuffer[y];
                    }
                }
            }

            //FFTs in z-direction
            AUTOPAS_OPENMP(parallel for schedule(static) firstprivate(inBuffer, outBuffer))
            for(unsigned y = 0; y < grid_dims[1]; y++){
                for(unsigned x = 0; x < grid_dims[0]; x++){
                    for(unsigned z = 0; z < grid_dims[2]; z++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        inBuffer[z] = fftBuffer3D[index1d];
                    }
                    forward(inBuffer, grid_dims[2], outBuffer);
                    for(unsigned z = 0; z < grid_dims[2]; z++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        out[index1d] = outBuffer[z];
                    }
                }
            }
        }

        /*void backward3DdirectionAt(ComplexGridType &in, ComplexGridType &out, std::array<unsigned int, 3> &grid_dims, int direction, int xIndx, int yIndx, int zIndx){
            for(unsigned int m = 0; m < grid_dims[direction]; m++){
                // copy values into buffer
                switch(direction){
                    case 0:
                        inBuffer[m] = in[m][yIndx][zIndx];
                        break;
                    case 1:
                        inBuffer[m] = in[xIndx][m][zIndx];
                        break;
                    case 2:
                        inBuffer[m] = in[xIndx][yIndx][m];
                        break;
                }
            }
            backward(inBuffer, grid_dims[direction], outBuffer);
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
            AUTOPAS_OPENMP(parallel for schedule(static) firstprivate(inBuffer, outBuffer))
            for (unsigned int z = 0; z < grid_dims[2]; z++){
                // do all FFTs for one x-Plane
                for(unsigned int y = 0; y < grid_dims[1]; y++){
                    // FFTs in x-direction
                    for(unsigned int x = 0; x < grid_dims[0]; x++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        inBuffer[x] = in[index1d];
                    }
                    backward(inBuffer, grid_dims[0], outBuffer);
                    for(unsigned int x = 0; x < grid_dims[0]; x++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        fftBuffer3D[index1d] = outBuffer[x];
                    }
                }
                for(unsigned int x = 0; x < grid_dims[0]; x++){
                    // FFTs in y-direction on the already transformed values
                    for(unsigned int y = 0; y < grid_dims[1]; y++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        inBuffer[y] = fftBuffer3D[index1d];
                    }
                    backward(inBuffer, grid_dims[1], outBuffer);
                    for(unsigned int y = 0; y < grid_dims[1]; y++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        fftBuffer3D[index1d] = outBuffer[y];
                    }
                }
            }

            //FFTs in z-direction
            AUTOPAS_OPENMP(parallel for schedule(static) firstprivate(inBuffer, outBuffer))
            for(unsigned y = 0; y < grid_dims[1]; y++){
                for(unsigned x = 0; x < grid_dims[0]; x++){
                    for(unsigned z = 0; z < grid_dims[2]; z++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        inBuffer[z] = fftBuffer3D[index1d];
                    }
                    backward(inBuffer, grid_dims[2], outBuffer);
                    for(unsigned z = 0; z < grid_dims[2]; z++){
                        unsigned int index1d = utils::ThreeDimensionalMapping::threeToOneD(x, y, z, grid_dims);
                        fftBuffer3D[index1d] = outBuffer[z];
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