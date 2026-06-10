#include <iostream>
#include <cmath>

// default size in x-direction
#define LX_DEFAULT 40.0
// default size in y-direction
#define LY_DEFAULT 20.0
// default number of particles in z-direction
// size of system in z = 2 + nz
#define NZ_DEFAULT 5

int main(int argv, char** argc)
{
    // if applicable read in size in x
    const double Lx = (argv >= 2)?atof(argc[1]):LX_DEFAULT;
    // if applicable read in size in y
    const double Ly = (argv >= 3)?atof(argc[2]):LY_DEFAULT;
    // if applicable read in number of particles in z
    const int nz = (argv >= 4)?atoi(argc[3]):NZ_DEFAULT;
    const double Lz = 2 + (double)nz;

    // output format: ASCII
    // <id> <x> <y> <z> <weight (=1.0 -> all values have the same weight)

    double l_a = std::sqrt(0.1) * Lx;
    double l_b = std::sqrt(0.1) * Lx;
    double l_c = 0.1 * Lx;
    double l_d = 0.3 * Lx;
    double l_e = 0.2 * Lx;
    double l_f = 0.3 * Lx;
    double l_g = 0.2 * Lx;

    int n_a = std::ceil(l_a); 
    int n_b = std::ceil(l_b);
    int n_c = 0.1 * Lx;
    int n_d = 0.3 * Lx;
    int n_e = 0.2 * Lx;
    int n_f = 0.3 * Lx;
    int n_g = 0.2 * Lx;

    double d_a = l_a / (double)n_a;
    double d_b = l_b / (double)n_b;
    double d_c = 1.0;
    double d_d = 1.0;
    double d_e = 1.0;
    double d_f = 1.0;
    double d_g = 1.0;

    int idx = 1;

    // create left part of 'A'
    for (int z = 1; z <= nz + 1; ++z)
    for (int t = 0; t <= n_a; ++t)
        std::cout << (idx++) << " "
                  << (1.0 +       (double)t/(double)n_a) * 0.1 * Lx << " "
                  << (1.0 + 3.0 * (double)t/(double)n_a) * 0.1 * Lx << " "
                  << (double)z << " " << 1.0 << std::endl;

    // create right part of 'A'
    for (int z = 1; z <= nz + 1; ++z)
    for (int t = 1; t <= n_b; ++t)
        std::cout << (idx++) << " "
                  << (2.0 +       (double)t/(double)n_b) * 0.1 * Lx << " "
                  << (4.0 - 3.0 * (double)t/(double)n_b) * 0.1 * Lx << " "
                  << (double)z << " " << 1.0 << std::endl;

    // create horizontal line of 'A'
    for (int z = 1; z <= nz + 1; ++z)
    for (int t = 1; t < n_c; ++t)
        std::cout << (idx++) << " "
                  << (1.5 +       (double)t/(double)n_c) * 0.1 * Lx << " "
                  << 2.5                                * 0.1 * Lx << " "
                  << (double)z << " " << 1.0 << std::endl;

    // create vertical part of first 'L'
    for (int z = 1; z <= nz + 1; ++z)
    for (int t = 1; t <= n_d; ++t)
        std::cout << (idx++) << " "
                  << 4.0                                * 0.1 * Lx << " "
                  << (1.0 + 3.0 * (double)t/(double)n_d) * 0.1 * Lx << " "
                  << (double)z << " " << 1.0 << std::endl;

    // create horizontal part of first 'L'
    for (int z = 1; z <= nz + 1; ++z)
    for (int t = 0; t <= n_e; ++t)
        std::cout << (idx++) << " "
                  << (4.0 + 2.0 * (double)t/(double)n_e) * 0.1 * Lx << " "
                  << 1.0                                * 0.1 * Lx << " "
                  << (double)z << " " << 1.0 << std::endl;

    // create vertical part of second 'L'
    for (int z = 1; z <= nz + 1; ++z)
    for (int t = 1; t <= n_f; ++t)
        std::cout << (idx++) << " "
                  << 7.0                                * 0.1 * Lx << " "
                  << (1.0 + 3.0 * (double)t/(double)n_f) * 0.1 * Lx << " "
                  << (double)z << " " << 1.0 << std::endl;

    // create horizontal part of second 'L'
    for (int z = 1; z <= nz + 1; ++z)
    for (int t = 0; t <= n_g; ++t)
        std::cout << (idx++) << " "
                  << (7.0 + 2.0 * (double)t/(double)n_g) * 0.1 * Lx << " "
                  << 1.0                                * 0.1 * Lx << " "
                  << (double)z << " " << 1.0 << std::endl;

}
