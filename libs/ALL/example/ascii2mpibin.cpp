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
    int dim[3];

    if (argc < 6)
    {
        std::cout << "usage: " << argv[0] << " <in_file (ASCII)> <out_file (binary)> <n_x> <n_y> <n_z>" << std::endl;
        MPI_Finalize();
        exit(-1);
    }

    // read dimension of process grid from command line
    dim[0] = atoi(argv[3]);
    dim[1] = atoi(argv[4]);
    dim[2] = atoi(argv[5]);

    int err;

    MPI_File outfile;

    std::ifstream infile;
    infile.open(argv[1]);
    
    err = MPI_File_open(MPI_COMM_WORLD, argv[2], MPI_MODE_CREATE | MPI_MODE_RDWR, MPI_INFO_NULL, &outfile);

    char line[256];

    long blockID;
    double values[4];
    int counter = 0;

    double max[3];

    // get maximum system size
    while(infile.getline(line,256))
    {
        if(infile.eof()) break;
        std::string str(line);
        std::istringstream istr(str);
        //std::cout << str << std::endl;
        istr >> blockID >> values[0] >> values[1] >> values[2] >> values[3];
        for (int i = 0; i < 3; ++i)
            if (max[i] < values[i]) max[i] = values[i];
    }        

    // go back to start of file
    infile.clear();
    infile.seekg(0,infile.beg);        
        
    double cell_size[3];

    for (int i = 0; i < 3; ++i)
    {
        // adjust the system size to the next multiple of 10
        max[i] = (double)((int)max[i] - ((int)max[i] % 10) + 10);
        cell_size[i] = max[i] / (double)dim[i];
    }

    std::cout << "domain layout: " << dim[0] << " " << dim[1] << " " << dim[2] << std::endl;
    std::cout << "system size: " << max[0] << " " << max[1] << " " << max[2] << std::endl;
    std::cout << "domain size: " << cell_size[0] << " " << cell_size[1] << " " << cell_size[2] << std::endl;

    int offset = 0;
    const int n_procs = dim[0] * dim[1] * dim[2];

    std::vector<long> loc_parts(dim[0]*dim[1]*dim[2]);
    std::vector<long> offset_domain(dim[0]*dim[1]*dim[2]);
    std::vector<int> curr_index(dim[0]*dim[1]*dim[2]);

    for (int i = 0; i < dim[0] * dim[1] * dim[2]; ++i)
    {
        loc_parts.at(i) = 0;
        offset_domain.at(i) = 0;
        curr_index.at(i) = 0;
    }

    // count particles per domain
    while(infile.getline(line,256))
    {
        if(infile.eof()) break;
        std::string str(line);
        std::istringstream istr(str);
        istr >> blockID >> values[0] >> values[1] >> values[2] >> values[3];

        int cell[3];

        for (int i = 0; i < 3; ++i)
            cell[i] = (int)(values[i] / cell_size[i]);

        loc_parts.at(cell[0]*dim[1]*dim[2]+cell[1]*dim[2]+cell[2])++;
    }

    std::partial_sum(loc_parts.begin(),loc_parts.end(),offset_domain.begin());

    // subtract the particles on own domain (partial sum includes local contribution)
    for (int i = 0; i < dim[0] * dim[1] * dim[2]; ++i)
    {
        offset_domain.at(i) -= loc_parts.at(i);
        int n_loc = (int)loc_parts.at(i);
        MPI_File_write_at_all(outfile,(MPI_Offset)(i*sizeof(long)),&offset_domain.at(i),1,MPI_LONG,MPI_STATUS_IGNORE);
        MPI_File_write_at_all(outfile,(MPI_Offset)(n_procs *sizeof(long) + i * sizeof(int)),&n_loc,1,MPI_INT,MPI_STATUS_IGNORE);
    }

    // go back to start of file
    infile.clear();
    infile.seekg(0,infile.beg);

    // copy point data to MPI I/O file
    while(infile.getline(line,256))
    {
        if(infile.eof()) break;
        std::string str(line);
        std::istringstream istr(str);
        istr >> blockID >> values[0] >> values[1] >> values[2] >> values[3];

        int cell[3];

        for (int i = 0; i < 3; ++i)
            cell[i] = (int)(values[i] / cell_size[i]);

        int domain = cell[0] * dim[1] * dim[2] + cell[1] * dim[2] + cell[2];

        long block_size = sizeof(long) + 4*sizeof(double);

        MPI_File_write_at_all(
                outfile,
                (MPI_Offset)
                   (n_procs*(sizeof(int) + sizeof(long))+(offset_domain.at(domain)+curr_index.at(domain))*block_size),
                &blockID,
                1,
                MPI_LONG,
                MPI_STATUS_IGNORE
        );

        MPI_File_write_at_all(
                outfile,
                (MPI_Offset)
                   (n_procs*(sizeof(int) + sizeof(long))+(offset_domain.at(domain)+curr_index.at(domain))*block_size+sizeof(long)),
                values,
                4,
                MPI_DOUBLE,
                MPI_STATUS_IGNORE
        );

        curr_index.at(domain)++;
    } 
    /*
    for (int ix = 0; ix < dim[0]; ++ix)
        for (int iy = 0; iy < dim[1]; ++iy)
            for (int iz = 0; iz < dim[2]; ++iz)
            {
                int index = iz + iy * dim[2] + ix * dim[1] * dim[2];
                int loc_count = 0;
                while(infile.getline(line,256))
                {
                    if(infile.eof()) break;
                    std::string str(line);
                    std::istringstream istr(str);
                    istr >> values[0] >> values[1] >> values[2] >> values[3] >> values[4];
                    
                    // if the point is within the current domain
                    if ( 
                            ((double)ix * cell_size[0]) <= values[1] && ((double)(ix + 1) * cell_size[0]) > values[1] &&
                            ((double)iy * cell_size[1]) <= values[2] && ((double)(iy + 1) * cell_size[1]) > values[2] &&
                            ((double)iz * cell_size[2]) <= values[3] && ((double)(iz + 1) * cell_size[2]) > values[3]
                       )
                    {
                        MPI_Offset write_location= (MPI_Offset)(5 * (offset + loc_count) * sizeof(double) + 2 * n_procs * sizeof(int));
                        MPI_File_write_at_all(outfile,write_location,values,5,MPI_DOUBLE,MPI_STATUS_IGNORE);
                        loc_count++;
                    }
                }
                std::cout << "index " << index << ": " << offset << " | " << loc_count << std::endl;
                MPI_File_write_at_all(outfile,(MPI_Offset)(index*sizeof(int)),&offset,1,MPI_INT,MPI_STATUS_IGNORE);
                MPI_File_write_at_all(outfile,(MPI_Offset)((n_procs + index)*sizeof(int)),&loc_count,1,MPI_INT,MPI_STATUS_IGNORE);
                offset += loc_count;
                infile.clear();
                infile.seekg(0,infile.beg);
            }
    */
    infile.close();
    MPI_File_close(&outfile);
    
    MPI_Finalize();
}
