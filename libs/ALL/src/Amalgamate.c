#include <stdlib.h>

#define STB_INCLUDE_IMPLEMENTATION
#define STB_INCLUDE_LINE_NONE
#include "stb_include.h"

int main(int argc, char **argv)
{
	if(argc<3)
	{
		printf("Usage: %s MAIN_FILE INCLUDE_PATH\n", argv[0]);
		return EXIT_FAILURE;
	}
	char Error[256];
	char *Amalgamated = stb_include_file(argv[1], "", argv[2], Error);
	if(!Amalgamated)
	{
		printf("%s", Error);
		return EXIT_FAILURE;
	}
	printf("%s",Amalgamated);
	return EXIT_SUCCESS;
}
