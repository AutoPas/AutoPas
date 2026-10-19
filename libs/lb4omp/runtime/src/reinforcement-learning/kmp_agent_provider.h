#ifndef AAA_H
#define AAA_H

#pragma once

#include <string>
#include <unordered_map>
#include <fstream>

#include "agents/kmp_agent.h"

class AgentProvider {
public:
    AgentProvider(const AgentProvider&) = delete;   // Disable copy constructor

    static AgentProvider& Get();

    /*
    * Creates and returns a pointer to a learning agent with the specified type and options.
    * */
    static Agent* create_agent(int type, std::string loopname, int states, int actions);

    /*
     * Main method for Reinforcement Learning AUTO method in LB4OMP.
     * */
    static int search(const std::string &loop_id, int agent_type, LoopData* stats, int dimension);

private:
    AgentProvider();
    ~AgentProvider();

    std::unordered_map<std::string, int> timesteps; // Save timestep progress
    std::unordered_map<std::string, Agent*> agents; // Save agent references across timesteps
    std::fstream ofs;                               // Single filestream for all the agents to write to

    static int* chunk_sizes;  // Stores the values of the chunk-sizes for the agent to try
    static std::fstream& get_filestream();
    static std::unordered_map<std::string, int>& get_timesteps();
    static std::unordered_map<std::string, Agent*>& get_agents();

    /*
     * Calculates the chunk sizes to try according to the formula:
     * chunk = iterations / (2^i * threads), where i also denotes the position
     * of the result in the array.
     * */
    static int calculate_chunks(int *array, int size, int n, int p);

    static void print_agent_params(Agent* agent, int* chunk_sizes_array, int dim);

    /*
     * Prints the stats of the agent to evaluate performance.
     * */
    static void print_agent_stats(Agent* agent, int cMethod);
};

#endif