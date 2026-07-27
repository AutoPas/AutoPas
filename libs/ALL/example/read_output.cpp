#include <stdlib.h>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>
#include <mpi.h>

int main (int argc, char** argv)
{
    MPI_Init(&argc, &argv);

    if (argc < 3)
    {
        std::cout << "usage: " << argv[0] << " <in_file (ASCII)> <n_p>" << std::endl;
        MPI_Finalize();
        exit(-1);
    }

    int n = atoi(argv[2]);

    MPI_File infile;
    int err;

    err = MPI_File_open(
                            MPI_COMM_WORLD,
                            argv[1],
                            MPI_MODE_CREATE | MPI_MODE_RDWR,
                            MPI_INFO_NULL,
                            &infile
                       );

    double positions[3];
    long blockID;
    long offset;
    int n_part;

    int blocksize = sizeof(long); // + 3 * sizeof(double);

    for (int d = 0; d < n; ++d)
    {
        MPI_File_read_at(
                            infile,
                            (MPI_Offset)(d * sizeof(int)),
                            &n_part,
                            1,
                            MPI_INT,
                            MPI_STATUS_IGNORE
                        );
        MPI_File_read_at(
                            infile,
                            (MPI_Offset)(n * sizeof(int) + d * sizeof(long)),
                            &offset,
                            1,
                            MPI_LONG,
                            MPI_STATUS_IGNORE
                        );
        for (int p = 0; p < n_part; ++p)
        {
            MPI_File_read_at(
                                infile,
                                (MPI_Offset)(
                                             offset + 
                                             p * blocksize 
                                            ),
                                &blockID,
                                1,
                                MPI_LONG,
                                MPI_STATUS_IGNORE
                            );
            /*
            MPI_File_read_at(
                                infile,
                                (MPI_Offset)(
                                             offset + 
                                             p * blocksize + 
                                             sizeof(long)
                                            ),
                                positions,
                                3,
                                MPI_DOUBLE,
                                MPI_STATUS_IGNORE
                            );
            */
            std::cout << "Rank: " << d
                      << " Offset: " << offset
                      << " N: " << n_part
                      << " Particle: " << p
                      << " ID: " << blockID
                      /*
                      << " Position: " << positions[0] << " "
                                       << positions[1] << " "
                                       << positions[2] << " "
                      */
                      << std::endl;
        }
    }

    MPI_Finalize();
}
