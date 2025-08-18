#pragma once

#include <complex>
#include <vector>
#include <math.h>
#include <array>

namespace autopas {
    
    class FFT {

        public:
        using RealGridType = typename std::vector<std::vector<std::vector<double>>>;
        using ComplexGridType = std::vector<std::vector<std::vector<std::complex<double>>>>;

        private:
        std::vector<std::complex<double>> inBuffer;
        std::vector<std::complex<double>> outBuffer;
        ComplexGridType fftBuffer3D;
        bool correct_scaling;

        void fft(std::vector<std::complex<double>> &x, unsigned int N, int start, std::vector<std::complex<double>> &result, bool forward, int stride = 1){
            if(N == 1){
                result[start] = x[start];
                return;
            }
            //split
            for(unsigned int i = 0; i < N/2; i++){
                result[start + i] = x[start + 2*i];
                result[start + N/2 + i] = x[start + 2*i+1];
            }
            fft(x, N/2, start, result, forward, 2*stride);
            fft(x, N/2, start + stride, result, forward, 2*stride);

            std::complex<double> w_n;
            if(forward){
                w_n = std::exp(std::complex(0.,-2*M_PI/ N));
            }else{
                w_n = std::exp(std::complex(0.,2*M_PI/ N));
            }
            std::complex<double> w = std::complex(1.,0.);

            for(unsigned int k = 0; k < N/2; k += stride){
                result[start + k] = result[start + k] + w * result[start + N/2 + k];
                result[start + (N/2) + k] = result[start + k] - w * result[start + N/2 + k];
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

            fftBuffer3D = ComplexGridType(N[0]);
            for(unsigned int i = 0; i < N[0]; i++){
                fftBuffer3D[i] = std::vector<std::vector<std::complex<double>>>(N[1]);
                for(unsigned int k = 0; k < N[1]; k++){
                    fftBuffer3D[i][k] = std::vector<std::complex<double>>(N[2]);
                }
            }
        }

        //TOOD move construction

        void forward(std::vector<std::complex<double>> &x, unsigned int N, std::vector<std::complex<double>> &result){
            fft(x, N, 0, result, true);
            if(correct_scaling){
                for(unsigned int i = 0; i < N; i++){
                    result[i] /= std::sqrt(2*M_PI);
                }
            }
        }

        void backward(std::vector<std::complex<double>> &x, unsigned int N, std::vector<std::complex<double>> &result){
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

        void forward3DdirectionAt(ComplexGridType &in, ComplexGridType &out, std::array<unsigned int, 3> &grid_dims, int direction, int xIndx, int yIndx, int zIndx){
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
        }

        void forward3D(RealGridType &in, ComplexGridType &out, std::array<unsigned int, 3> &grid_dims){
            AUTOPAS_OPENMP(parallel for schedule(static))
            for(unsigned int i = 0; i < grid_dims[0]; i++){
                for(unsigned int j = 0; j < grid_dims[1]; j++){
                    for(unsigned int k = 0; k < grid_dims[2]; k++){
                        fftBuffer3D[i][j][k] = std::complex<double>(in[i][j][k]);
                    }
                }
            }

            AUTOPAS_OPENMP(parallel for schedule(static) firstprivate(inBuffer, outBuffer))
            for (unsigned int i = 0; i < grid_dims[0]; i++){
                // do all FFTs for one x-Plane
                for(unsigned int j = 0; j < grid_dims[1]; j++){
                    // FFTs in z-direction
                    for(unsigned int k = 0; k < grid_dims[2]; k++){
                        inBuffer[k] = fftBuffer3D[i][j][k];
                    }
                    forward(inBuffer, grid_dims[2], outBuffer);
                    for(unsigned int k = 0; k < grid_dims[2]; k++){
                        fftBuffer3D[i][j][k] = outBuffer[k];
                    }
                }
                for(unsigned int k = 0; k < grid_dims[2]; k++){
                    // FFTs in y-direction on the already transformed values
                    for(unsigned int j = 0; j < grid_dims[1]; j++){
                        inBuffer[j] = fftBuffer3D[i][j][k];
                    }
                    forward(inBuffer, grid_dims[1], outBuffer);
                    for(unsigned int j = 0; j < grid_dims[1]; j++){
                        fftBuffer3D[i][j][k] = outBuffer[j];
                    }
                }
            }

            //FFTs in x-direction
            AUTOPAS_OPENMP(parallel for schedule(static) firstprivate(inBuffer, outBuffer))
            for(unsigned j = 0; j < grid_dims[1]; j++){
                for(unsigned k = 0; k < grid_dims[2]; k++){
                    for(unsigned i = 0; i < grid_dims[0]; i++){
                        inBuffer[i] = fftBuffer3D[i][j][k];
                    }
                    forward(inBuffer, grid_dims[0], outBuffer);
                    for(unsigned i = 0; i < grid_dims[0]; i++){
                        out[i][j][k] = outBuffer[i];
                    }
                }
            }
        }

        void backward3DdirectionAt(ComplexGridType &in, ComplexGridType &out, std::array<unsigned int, 3> &grid_dims, int direction, int xIndx, int yIndx, int zIndx){
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
        }

        void backward3D(ComplexGridType &in, RealGridType &out, std::array<unsigned int, 3> &grid_dims){
            AUTOPAS_OPENMP(parallel for schedule(static) firstprivate(inBuffer, outBuffer))
            for (unsigned int i = 0; i < grid_dims[0]; i++){
                // do all FFTs for one x-Plane
                for(unsigned int j = 0; j < grid_dims[1]; j++){
                    // FFTs in z-direction
                    for(unsigned int k = 0; k < grid_dims[2]; k++){
                        inBuffer[k] = in[i][j][k];
                    }
                    backward(inBuffer, grid_dims[2], outBuffer);
                    for(unsigned int k = 0; k < grid_dims[2]; k++){
                        fftBuffer3D[i][j][k] = outBuffer[k];
                    }
                }
                for(unsigned int k = 0; k < grid_dims[2]; k++){
                    // FFTs in y-direction on the already transformed values
                    for(unsigned int j = 0; j < grid_dims[1]; j++){
                        inBuffer[j] = fftBuffer3D[i][j][k];
                    }
                    backward(inBuffer, grid_dims[1], outBuffer);
                    for(unsigned int j = 0; j < grid_dims[1]; j++){
                        fftBuffer3D[i][j][k] = outBuffer[j];
                    }
                }
            }

            //FFTs in x-direction
            AUTOPAS_OPENMP(parallel for schedule(static) firstprivate(inBuffer, outBuffer))
            for(unsigned j = 0; j < grid_dims[1]; j++){
                for(unsigned k = 0; k < grid_dims[2]; k++){
                    for(unsigned i = 0; i < grid_dims[0]; i++){
                        inBuffer[i] = fftBuffer3D[i][j][k];
                    }
                    backward(inBuffer, grid_dims[0], outBuffer);
                    for(unsigned i = 0; i < grid_dims[0]; i++){
                        fftBuffer3D[i][j][k] = outBuffer[i];
                    }
                }
            }

            // Make output real
            AUTOPAS_OPENMP(parallel for schedule(static))
            for(unsigned int i = 0; i < grid_dims[0]; i++){
                for(unsigned int j = 0; j < grid_dims[1]; j++){
                    for(unsigned int k = 0; k < grid_dims[2]; k++){
                        out[i][j][k] = std::real(fftBuffer3D[i][j][k]);
                    }
                }
            }
        }

        
    };
}