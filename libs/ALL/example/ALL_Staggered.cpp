/*
   Copyright 2020-2020 Stephan Schulz, Forschungszentrum Juelich GmbH, Germany
   Copyright 2018-2024 Rene Halver, Forschungszentrum Juelich GmbH, Germany
   Copyright 2018-2020 Godehard Sutmann, Forschungszentrum Juelich GmbH, Germany

   Redistribution and use in source and binary forms, with or without
   modification, are permitted provided that the following conditions are met:

   1. Redistributions of source code must retain the above copyright notice,
   this list of conditions and the following disclaimer.

   2. Redistributions in binary form must reproduce the above copyright notice,
   this list of conditions and the following disclaimer in the documentation
   and/or other materials provided with the distribution.

   3. Neither the name of the copyright holder nor the names of its contributors
   may be used to endorse or promote products derived from this software without
   specific prior written permission.

   THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
   AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
   IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
   ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
   LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
   CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
   SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
   INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
   CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
   ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
   POSSIBILITY OF SUCH DAMAGE.
 */

#include <stdlib.h>
#include <stdio.h>
// #include <time.h>
#include <mpi.h>
#include <string>
#include <fstream>
#include <iostream>

// #define ALL_VTK_OUTPUT

#ifdef USE_AMALGAMATED
#include "ALL_Amalgamated.hpp"
#else
#include <ALL.hpp>
#endif

// Run fun in order of ranks
// Todo(s.schulz): This seems to only work roughly with the result width an 32 ranks, with up to 16 it seems to work correctly.
//                 Adding sleep(1) also orders everything correctly. So this is probably a flushing problem.
//                 It also exists for the cout stream with endl.
#define MPI_RUN_ORDER(comm, rank, max_ranks, fun)               \
	{                                                           \
		int MPI_RO_IT;                                          \
		for (MPI_RO_IT = 0; MPI_RO_IT < max_ranks; MPI_RO_IT++) \
		{                                                       \
			if (MPI_RO_IT == rank)                              \
			{                                                   \
				fun;                                            \
				MPI_Barrier(comm);                              \
			}                                                   \
			else                                                \
			{                                                   \
				MPI_Barrier(comm);                              \
			}                                                   \
		}                                                       \
	}

// Quick and dirty helper function. Assumes comm, rank and max_ranks
// CAVEAT: the function call must be wrapped in () if it contains a comma
#define MPI_RUN_ORDER_DEF(fun) MPI_RUN_ORDER(MPI_COMM_WORLD, MyRank, MaximumRank, fun)

// void millisleep(unsigned long ms)
//{
//	struct timespec SleepTime;
//	SleepTime.tv_sec = ms/1000;
//	SleepTime.tv_nsec = (ms%1000)*1000000L;
//	struct timespec RemainingTime;
//	while(nanosleep(&SleepTime, &RemainingTime))
//	{
//		SleepTime.tv_sec = RemainingTime.tv_sec;
//		SleepTime.tv_nsec = RemainingTime.tv_nsec;
//	}
// }

void print_width(int rank, double width, double bottom, double top)
{
	printf("[%03d] Result Width: %10.6f (%10.6f -- %10.6f)\n", rank, width, bottom, top);
	fflush(stdout);
}

void print_testing_output(int rank, std::vector<ALL::Point<double>> &vertices, int timestep)
{
	//	printf("[%4d,%03d] Result Width: %10.6f %10.6f %10.6f\n",
	//			timestep,
	//			rank,
	//			vertices.at(1)[0]-vertices.at(0)[0],
	//			vertices.at(1)[1]-vertices.at(0)[1],
	//			vertices.at(1)[2]-vertices.at(0)[2]);
	//	fflush(stdout);
	printf("[%4d,%03d,%02d] Result Vertex: %10.6f %10.6f %10.6f\n",
		   timestep,
		   rank,
		   0,
		   vertices.at(0)[0],
		   vertices.at(0)[1],
		   vertices.at(0)[2]);
	fflush(stdout);
	printf("[%4d,%03d,%02d] Result Vertex: %10.6f %10.6f %10.6f\n",
		   timestep,
		   rank,
		   1,
		   vertices.at(1)[0],
		   vertices.at(1)[1],
		   vertices.at(1)[2]);
	fflush(stdout);
}

void print_loc(int rank, int *loc, int *size)
{
	printf("[%03d] Location: (%3d,%3d,%3d)/(%3d,%3d,%3d)\n", rank, loc[0], loc[1], loc[2], size[0], size[1], size[2]);
	fflush(stdout);
}

void print_domain(int rank, double *verts)
{
	printf("[%03d] Lower: %g\t%g\t%g\n", rank, verts[0], verts[1], verts[2]);
	printf("[%03d] Upper: %g\t%g\t%g\n", rank, verts[3], verts[4], verts[5]);
	fflush(stdout);
}

void print_work(int rank, double work)
{
	printf("[%03d] Work: %g\n", rank, work);
	fflush(stdout);
}

void print_binary(int step, int rank, double work, std::vector<ALL::Point<double>> &vertices, int *loc, int *size, MPI_File fh)
{
	int nranks = size[0] * size[1] * size[2];

	MPI_Offset offset = (step * size[0] * size[1] * size[2] * 11 + rank * 11) * sizeof(double);

	double buf[11];

	buf[0] = (double)rank;
	buf[1] = work;
	buf[2] = vertices.at(0)[0];
	buf[3] = vertices.at(0)[1];
	buf[4] = vertices.at(0)[2];
	buf[5] = vertices.at(1)[0];
	buf[6] = vertices.at(1)[1];
	buf[7] = vertices.at(1)[2];
	buf[8] = (double)(loc[0]);
	buf[9] = (double)(loc[1]);
	buf[10] = (double)(loc[2]);

	MPI_Status state;
	MPI_File_write_at(fh, offset, buf, 11, MPI_DOUBLE, &state);
}

void convert_verts(std::vector<ALL::Point<double>> *vv, double *verts)
{
	verts[0] = vv->at(0)[0];
	verts[1] = vv->at(0)[1];
	verts[2] = vv->at(0)[2];
	verts[3] = vv->at(1)[0];
	verts[4] = vv->at(1)[1];
	verts[5] = vv->at(1)[2];
}

int main(int argc, char **argv)
{
	MPI_Init(&argc, &argv);
	int CurrentStep = 0;
	const int NumberOfSteps = 50;

	// set the output directory for the binary file
	std::string outputDir;
	if (argc > 1)
	{
		outputDir.assign(argv[1]);
	}
	else
		outputDir.assign(".");

	const int Dimensions = 3;
	const int LoadbalancerGamma = 0; // ignored for staggered method
	// Todo(s.schulz): is gamma also ignored for tensor method?
	// Todo(s.schulz): maybe change to automatic destruction on scope end
	// The ALL::TENSOR method can be used just like the staggered method.
	// Just exchange ALL::STAGGERED with ALL::TENSOR in the constructor.
	ALL::ALL<double, double> *jall = new ALL::ALL<double, double>(ALL::STAGGERED, Dimensions, LoadbalancerGamma);
	int MyLocation[3] = {0};
	// All domains are placed along a line in z direction, even though they are three dimensional
	MPI_Comm_rank(MPI_COMM_WORLD, &MyLocation[2]);
	int MyRank = MyLocation[2];

	int NumberOfProcesses[3] = {1, 1, 1};
	MPI_Comm_size(MPI_COMM_WORLD, &NumberOfProcesses[2]);
	int MaximumRank = NumberOfProcesses[2];

	std::stringstream ss;
	ss << outputDir << "/" << "ALL_Staggered_" << NumberOfProcesses[2] << ".bin";
	std::string file = ss.str();
	MPI_File fh;
	MPI_File_delete(file.c_str(), MPI_INFO_NULL);
	MPI_File_open(MPI_COMM_WORLD, file.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh);
	MPI_File_set_size(fh, NumberOfSteps * NumberOfProcesses[2] * 11 * sizeof(MPI_DOUBLE));
	MPI_File_preallocate(fh, NumberOfSteps * NumberOfProcesses[2] * 11 * sizeof(MPI_DOUBLE));

	if (MyRank == 0)
	{
		printf("Ranks: %d\nNumber of Steps: %d\n", MaximumRank, NumberOfSteps);
		fflush(stdout);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_RUN_ORDER_DEF((print_loc(MyRank, MyLocation, NumberOfProcesses)));

	// For a cartesian communicator this is not required, but we are using
	// MPI_COMM_WORLD here.
	std::vector<int> MyLocationVector(MyLocation, MyLocation + 3);
	std::vector<int> NumberOfProcessesVector(NumberOfProcesses, NumberOfProcesses + 3);
	jall->setProcGridParams(MyLocationVector, NumberOfProcessesVector);

	// If we want to set a minimum domain size:
	std::vector<double> MinimumDomainSize{0.1, 0.1, 0.1};
	jall->setMinDomainSize(MinimumDomainSize);

	jall->setCommunicator(MPI_COMM_WORLD);

	// We also set the optional process tag for the output.
	// This can be useful if we want to know which of 'our'
	// ranks is which in the output produces by the library.
	// The ranks used inside the library do not necessarily
	// match our own numbering.
	jall->setProcTag(MyRank);

	jall->setup();
	// Todo(s.schulz): document: what exactly must be set before setup()?

	// A first domain distribution must be given to the balancer.
	// We use the provided ALL::Point class to define the vertices,
	// but a simple double array can also be used. We need 2 vertices
	// which correspond to the two opposing corners.
	std::vector<ALL::Point<double>> DomainVertices(2, ALL::Point<double>(3));
	const double DomainSize = 1.0; // Domain size
	// We create a cubic domain initially
	for (int VertexIndex = 0; VertexIndex < 2; VertexIndex++)
	{
		for (int DimensionIndex = 0; DimensionIndex < Dimensions; DimensionIndex++)
		{
			DomainVertices.at(VertexIndex)[DimensionIndex] = (MyLocation[DimensionIndex] + VertexIndex) * DomainSize;
		}
	}
	double VertexArray[6];
	convert_verts(&DomainVertices, VertexArray);
	MPI_RUN_ORDER_DEF((print_domain(MyRank, VertexArray)));
	jall->setVertices(DomainVertices);

	// Calculate the work of our domain. Here we just use
	double MyWork = (double)MyRank + 1.;
	jall->setWork(MyWork);
	MPI_RUN_ORDER_DEF((print_work(MyRank, MyWork)));
	for (CurrentStep = 0; CurrentStep < NumberOfSteps; CurrentStep++)
	{
		// In a real code we need to set the updated work in each
		// iteration before calling balance()
		if (MyRank == 0)
		{
			printf("Starting step: %d/%d\n", CurrentStep + 1, NumberOfSteps);
			fflush(stdout);
		}
#ifdef ALL_VTK_OUTPUT_EXAMPLE
		try
		{
			jall->printVTKoutlines(CurrentStep);
		}
		catch (ALL::FilesystemErrorException &e)
		{
			std::cout << e.what() << std::endl;
			std::exit(2);
			// Maybe also treat this error in some way, e.g. whether the simulation
			// should abort if no output could be written or not.
		}
#endif
		jall->balance();

		std::vector<ALL::Point<double>> NewVertices = jall->getVertices();
		// MPI_RUN_ORDER(MPI_COMM_WORLD, MyRank, MaximumRank, (print_width(MyRank, NewVertices.at(1)[2]-NewVertices.at(0)[2], NewVertices.at(0)[2], NewVertices.at(1)[2])));
		// MPI_RUN_ORDER_DEF((print_testing_output(MyRank, NewVertices, CurrentStep + 1)));
		print_binary(CurrentStep, MyRank, MyWork, NewVertices, MyLocation, NumberOfProcesses, fh);

		// Maybe print our new domain? Or not..
		// convert_verts(&NewVertices, VertexArray);
		// MPI_RUN_ORDER_DEF((print_domain(MyRank, VertexArray)));
		// jall->getWork(MyWork);
		// MPI_RUN_ORDER_DEF((print_work(MyRank,MyWork)));
		// MPI_Barrier(MPI_COMM_WORLD);
	}
#ifdef ALL_VTK_OUTPUT_EXAMPLE
	try
	{
		jall->printVTKoutlines(CurrentStep);
	}
	catch (ALL::FilesystemErrorException &e)
	{
		std::cout << e.what() << std::endl;
	}
#endif

	delete jall;
	MPI_File_close(&fh);

	/* Check correctness of binary file (if in doubt)*/
	/*
	if (MyRank == 0)
	{
		std::string filename = ss.str();
		std::ifstream infile(filename, std::ios::binary);

		double d;
		for (int nSteps = 0; nSteps < NumberOfSteps; ++nSteps)
			for (int n = 0; n < NumberOfProcesses[2]; ++n)
			{
				for (int i = 0; i < 11; ++i)
				{
					infile.read(reinterpret_cast<char *>(&d), sizeof d);
					std::cerr << d << " ";
				}
				std::cerr << '\n';
			}
		infile.close();
	}
	*/

	MPI_Finalize();
	return EXIT_SUCCESS;
}