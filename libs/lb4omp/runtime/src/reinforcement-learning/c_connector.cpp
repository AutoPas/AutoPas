#include <cstdlib>
#include <iostream>

#include "c_connector.h"
#include "kmp_agent_provider.h"
#include "kmp_loopdata.h"
#include "utils/utils.h"


#ifdef __cplusplus
extern "C" {
#endif

static AgentProvider *agent_instance = NULL;

int search(const char *loop_id, int agent_type, LoopData *stats, int dimension) {
    return AgentProvider::search(loop_id, agent_type, stats, dimension);
}

void setEnvVar(const char *var_name, int value) {
    set_env_var(var_name, value);
}

void setLoopData(const char *loop_id, LoopData *loop_data) {
    set_loop_data(loop_id, loop_data);
}

LoopData *getLoopData(const char *loop_id) {
    return get_loop_data(loop_id);
}

#ifdef __cplusplus
}
#endif