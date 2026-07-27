/*
   Copyright 2020-2020 Stephan Schulz, Forschungszentrum Juelich GmbH, Germany

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
#include <string.h>
#include <math.h>

#define STB_DEFINE
#define STB_DEBUG
#include "stb_arr.h"

// Returns 0 on Match and 2 on Bad Match.
// Other errors return different codes.

#define EXIT_NOMATCH 2

// This assumes the following format:
// First two lines .must. be:
//  Ranks: %d
//  Number of Steps: %d
// Then lines are ignored except
//  [%d,%d,%d] Result Vertex: %f %f %f
// With: [Step, Rank, Vertex]
// The lines need not be ordered, but must all be placed after Number
// of Steps
// Output: Data[Step][Rank][Vertex][Dimension]
double ****LoadTest(char *FileName)
{
	FILE* File = fopen(FileName,"r");
	if(!File)
	{
		fprintf(stderr,"Could not open file '%s'\n", FileName);
		exit(EXIT_FAILURE);
	}
	const int LineLength = 8192;
	char Line[LineLength];
	int Dimension = 3;
	if(fgets(Line, LineLength, File) == NULL)
	{
		fprintf(stderr,"Could not read file '%s'\n", FileName);
		exit(EXIT_FAILURE);
	}
	int Ranks;
	int ReadCount = sscanf(Line, "Ranks: %d", &Ranks);
	if(ReadCount<1)
	{
		fprintf(stderr,"Could not parse ranks in '%s'\n", FileName);
		exit(EXIT_FAILURE);
	}
	if(fgets(Line, LineLength, File) == NULL)
	{
		fprintf(stderr,"Could not read file '%s'\n", FileName);
		exit(EXIT_FAILURE);
	}
	int Steps;
	ReadCount = sscanf(Line, "Number of Steps: %d", &Steps);
	if(ReadCount<1)
	{
		fprintf(stderr,"Could not parse number of steps in '%s'\n", FileName);
		exit(EXIT_FAILURE);
	}
	double ****Data = NULL;
	stb_arr_setlen(Data, Steps);
	memset(Data, 0, sizeof(Data)*Steps);
	int CurrentStep = 0;
	while(fgets(Line, LineLength, File) != NULL)
	{
		double InRankData[3];
		int CurrentRank;
		int CurrentVertex;
		ReadCount = sscanf(Line, "[%d,%d,%d] Result Vertex: %lf %lf %lf",
				&CurrentStep,
				&CurrentRank,
				&CurrentVertex,
				&InRankData[0], &InRankData[1], &InRankData[2]);
		if(ReadCount!=6) continue;
		CurrentStep--; //The input is 1 indexed, the code 0 indexed.
		double ***StepData = Data[CurrentStep];
		if(!stb_arr_valid(StepData, CurrentRank))
		{
			// increase array and initialise values
			int CurrentLength = stb_arr_len(StepData);
			stb_arr_setlen(StepData, CurrentRank+1);
			memset(&StepData[CurrentLength], 0, sizeof(StepData)*(CurrentRank+1-CurrentLength));
			Data[CurrentStep] = StepData;
		}
		double **RankData = StepData[CurrentRank];
		if(!stb_arr_valid(RankData, CurrentVertex))
		{
			int CurrentLength = stb_arr_len(RankData);
			stb_arr_setlen(RankData, CurrentVertex+1);
			memset(&RankData[CurrentVertex], 0, sizeof(RankData)*(CurrentVertex+1-CurrentLength));
			StepData[CurrentRank] = RankData;
		}
		double *VertexData = RankData[CurrentVertex];
		stb_arr_setlen(VertexData, Dimension);
		RankData[CurrentVertex] = VertexData;
		VertexData[0] = InRankData[0];
		VertexData[1] = InRankData[1];
		VertexData[2] = InRankData[2];
		//printf("IN: %s", Line);
		//printf("ME: [%4d,%03d,%02d] Result Vertex: %10.6f %10.6f %10.6f\n",
		//		CurrentStep+1,
		//		CurrentRank,
		//		CurrentVertex,
		//		InRankData[0], InRankData[1], InRankData[2]);
	}
	fclose(File);
	return Data;
}

// return 0 on difference and >1 if same
int CompareTests(double ****Test1, double ****Test2)
{
	double MaximumDifference = 0;
	int CurrentStep, CurrentRank, CurrentVertex, i;
	if(stb_arr_len(Test1) != stb_arr_len(Test2)){
		printf("Different number of steps: %d %d\n", stb_arr_len(Test1), stb_arr_len(Test2));
		return 0;
	}
	//printf("Steps: %d %d\n", stb_arr_len(Test1), stb_arr_len(Test2));
	for(CurrentStep=0; CurrentStep<stb_arr_len(Test1); CurrentStep++)
	{
		if(stb_arr_len(Test1[CurrentStep]) != stb_arr_len(Test2[CurrentStep]))
		{
			printf("Different number of ranks in Step %d: %d %d\n",
					CurrentStep,
					stb_arr_len(Test1[CurrentStep]),
					stb_arr_len(Test2[CurrentStep]));
			return 0;
		}
		//printf(" Step: %3d  Ranks: %d %d\n", CurrentStep, stb_arr_len(Test1[CurrentStep]), stb_arr_len(Test2[CurrentStep]));
		for(CurrentRank=0; CurrentRank<stb_arr_len(Test1[CurrentStep]); CurrentRank++)
		{
			if(stb_arr_len(Test1[CurrentStep][CurrentRank]) != stb_arr_len(Test2[CurrentStep][CurrentRank]))
			{
				printf("Different number of vertices in Step %d and Rank %d: %d %d\n",
						CurrentStep,
						CurrentRank,
						stb_arr_len(Test1[CurrentStep][CurrentRank]),
						stb_arr_len(Test2[CurrentStep][CurrentRank]));
				return 0;
			}
			for(CurrentVertex=0; CurrentVertex<stb_arr_len(Test1[CurrentStep][CurrentRank]); CurrentVertex++)
			{
				assert(stb_arr_valid(Test2[CurrentStep][CurrentRank], CurrentVertex));
				for(i=0; i<stb_arr_len(Test1[CurrentStep][CurrentRank][CurrentVertex]); i++)
				{
					double Value1 = Test1[CurrentStep][CurrentRank][CurrentVertex][i];
					double Value2 = Test2[CurrentStep][CurrentRank][CurrentVertex][i];
					double AbsVal = fabs(Value1-Value2);
					MaximumDifference = MaximumDifference>AbsVal?MaximumDifference:AbsVal;
					if(AbsVal>0.001) printf("Found deviation of %g at Step: %u Rank: %u Vertex: %u, between %g and %g\n",
							AbsVal,
							CurrentStep,
							CurrentRank,
							CurrentVertex,
							Value1, Value2);
				}
			}
		}
	}
	printf("MaximumDifference: %g\n", MaximumDifference);
	return MaximumDifference<0.001;
}

int main(int argc, char **argv)
{
	if(argc<3)
	{
		printf("Usage: %s KNOWN_GOOD TEST_OUTPUT\n", argv[0]);
		return EXIT_FAILURE;
	}
	char *GoodFileName = argv[1];
	char *TestFileName = argv[2];

	// TestData2[Step][Rank][Vertex][Dimension]
	double ****GoodData = LoadTest(GoodFileName);
	double ****TestData = LoadTest(TestFileName);

	int Matches = CompareTests(GoodData, TestData);

	if(Matches)
		return EXIT_SUCCESS;
	else
		return EXIT_NOMATCH;
}
