#include <iostream>

#include "fortran_connector.h"
#include "kmp_agent_provider.h"
#include "kmp_loopdata.h"

#ifdef __cplusplus
extern "C" {
#endif

int my_cpp_function(int arg) {
    std::cout << "Hello from C++! The argument is: " << arg << std::endl;
    arg = arg * 2;
    return arg;
}

static AgentProvider *agent_instance = NULL;

int search(const char *loop_id, int agent_type, LoopData *stats, int dimension) {
    return AgentProvider::search(loop_id, agent_type, stats, dimension);
}

#ifdef __cplusplus
}
#endif