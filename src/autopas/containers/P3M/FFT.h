#pragma once

#include <complex>
#include <vector>
#include <numbers>
#include <array>

namespace autopas {
    
    class FFT {

        using RealGridType = typename std::vector<std::vector<std::vector<double>>>;
        using ComplexGridType = std::vector<std::vector<std::vector<std::complex<double>>>>;

        double pi = 3.1415926535897932384626433832795028841971693993751058209749445923078164062;
        std::vector<std::complex<double>> inBuffer;
        std::vector<std::complex<double>> outBuffer;
        std::vector<std::complex<double>> fftBuffer;
        bool correct_scaling;

        private:
        void fft(std::vector<std::complex<double>> &x, int N, int start, std::vector<std::complex<double>> &result, bool forward){
            if(N == 1){
                result[start] = x[start];
                return;
            }
            //split
            for(int i = 0; i < N/2; i++){
                result[start + i] = x[start + 2*i];
                result[start + N/2 + i] = x[start + 2*i+1];
            }
            fft(result, N/2, start, x, forward);
            fft(result, N/2, start + N/2, x, forward);

            std::complex<double> w_n;
            if(forward){
                w_n = std::exp(std::complex(0.,2*pi/ N));
            }else{
                w_n = std::exp(std::complex(0.,-2*pi/ N));
            }
            std::complex<double> w = std::complex(1.,0.);

            for(int k = 0; k < N/2; k++){
                result[start + k] = x[start + k] + w * x[start + N/2 + k];
                result[start + (N/2) + k] = x[start + k] - w * x[start + N/2 + k];
                w *=  w_n;
            }
        }

        public:

        FFT(){}

        FFT(int bufferSize, bool correct_scaling = false) : inBuffer(std::vector<std::complex<double>>(bufferSize)), outBuffer(std::vector<std::complex<double>>(bufferSize)), fftBuffer(std::vector<std::complex<double>>(bufferSize)), correct_scaling(correct_scaling){}

        //TOOD move construction

        void forward(std::vector<std::complex<double>> &x, int N, std::vector<std::complex<double>> &result){
            fft(x, N, 0, result, true);
            if(correct_scaling){
                for(int i = 0; i < N; i++){
                    result[i] /= std::sqrt(2*pi);
                }
            }
        }

        void backward(std::vector<std::complex<double>> &x, int N, std::vector<double> &result){
            fft(x, N, 0, fftBuffer, false);
            if(correct_scaling){
                for(int i = 0; i < N; i++){
                    result[i] = std::real(fftBuffer[i]) * std::sqrt(2*pi) / N;
                }
            }else{
                for(int i = 0; i < N; i++){
                    result[i] = std::real(fftBuffer[i]) / N;
                }
            }
            
        }

        void forward3Ddirection(RealGridType &in, ComplexGridType &out, std::array<int, 3> &grid_dims, int direction){
            for(int k = 0; k < grid_dims[(direction + 1) % 3]; k++){
                for(int l = 0; l < grid_dims[(direction + 2) % 3]; l++){
                    // ffts in z-direction
                    for(int m = 0; m < grid_dims[direction]; m++){
                        // copy values into buffer
                        switch(direction){
                            case 0:
                                inBuffer[m] = std::complex<double>(in[m][k][l]);
                                break;
                            case 1:
                                inBuffer[m] = std::complex<double>(in[l][m][k]);
                                break;
                            case 2:
                                inBuffer[m] = std::complex<double>(in[k][l][m]);
                                break;
                        }
                        
                    }
                    forward(inBuffer, grid_dims[direction], outBuffer);
                    for(int m = 0; m < grid_dims[2]; m++){
                        // copy values into buffer
                        switch(direction){
                            case 0:
                                out[m][k][l] += outBuffer[m];
                                break;
                            case 1:
                                out[l][m][k] += outBuffer[m];
                                break;
                            case 2:
                                out[k][l][m] += outBuffer[m];
                                break;
                        }
                    }
                }
            }
        }

        void forward3D(RealGridType &in, ComplexGridType &out, std::array<int, 3> &grid_dims){
            for(int i = 0; i < 3; i++){
                forward3Ddirection(in, out, grid_dims, i);
            }

            /*for(int xi = 0; xi < grid_dims[0]; xi++){
                for(int yi = 0; yi < grid_dims[1]; yi++){
                    // ffts in z-direction
                    for(int zi = 0; zi < grid_dims[2]; zi++){
                        // copy values into buffer over
                        realBuffer[zi] = in[xi][yi][zi];
                    }
                    forward(realBuffer, grid_dims[2], complexBuffer);
                    for(int zi = 0; zi < grid_dims[2]; zi++){
                        // copy values into buffer over
                        out[xi][yi][zi] += complexBuffer[zi];
                    }
                }
            }*/
        }

        void backward3Ddirection(ComplexGridType &in, RealGridType &out, std::array<int, 3> &grid_dims, int direction){
            for(int k = 0; k < grid_dims[(direction + 1) % 3]; k++){
                for(int l = 0; l < grid_dims[(direction + 2) % 3]; l++){
                    // ffts in z-direction
                    for(int m = 0; m < grid_dims[direction]; m++){
                        // copy values into buffer
                        switch(direction){
                            case 0:
                                inBuffer[m] = in[m][k][l];
                                break;
                            case 1:
                                inBuffer[m] = in[l][m][k];
                                break;
                            case 2:
                                inBuffer[m] = in[k][l][m];
                                break;
                        }
                        
                    }
                    forward(inBuffer, grid_dims[direction], outBuffer);
                    for(int m = 0; m < grid_dims[2]; m++){
                        // copy values into buffer
                        switch(direction){
                            case 0:
                                out[m][k][l] += std::real(outBuffer[m]);
                                break;
                            case 1:
                                out[l][m][k] += std::real(outBuffer[m]);
                                break;
                            case 2:
                                out[k][l][m] += std::real(outBuffer[m]);
                                break;
                        }
                    }
                }
            }
        }

        void backward3D(ComplexGridType &in, RealGridType &out, std::array<int, 3> &grid_dims){
            for(int i = 0; i < 3; i++){
                backward3Ddirection(in, out, grid_dims, i);
            }
        }

        
    };
}