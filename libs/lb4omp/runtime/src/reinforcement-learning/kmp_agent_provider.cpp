// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//

#include <string>
#include <iostream>
#include <iomanip>

#include "utils/utils.h"
#include "kmp_agent_provider.h"
#include "initializers/kmp_init_type.h"
#include "policies/kmp_policy_type.h"
#include "agents/kmp_agent.h"
#include "agents/kmp_q_learner.h"
#include "agents/kmp_sarsa_learner.h"
#include "agents/kmp_expected_sarsa_learner.h"
#include "agents/kmp_doubleq_learner.h"
#include "agents/kmp_qv_learner.h"
#include "agents/kmp_r_learner.h"
#include "agents/kmp_deepq_learner.h"
#include "agents/defaults.h"


int *AgentProvider::chunk_sizes;  // Stores the values of the chunk-sizes for the agent to try

AgentProvider &AgentProvider::Get() {
    static AgentProvider instance;
    return instance;
}

AgentProvider::AgentProvider() {
    std::string fileData{defaults::KMP_RL_AGENT_STATS};
    read_env_string("KMP_RL_AGENT_STATS", fileData);
    if (fileData != "notOutputtingStats") {
        ofs.open(fileData, std::ofstream::out | std::ofstream::app);
        ofs << "loop_id,\tstep,\tmethod,\talpha,\t\teps,\t\treward,\t\tlow,\t\thigh,\t\tcurrent" << std::endl;
    }
}

AgentProvider::~AgentProvider() {
    ofs.close();
}

// public
int AgentProvider::search(const std::string &loop_id, int agent_type, LoopData *stats, int dimension) {
#if (_RL_DEBUG > 0)
    std::cout << "[AgentProvider::search] Loop: " << loop_id << std::endl;
    std::cout << "[AgentProvider::search] Current LoopTime: " << stats->cTime << std::endl;
#endif
    if (!AgentProvider::get_timesteps().count(loop_id)) {
        // Enters if the agent does not exist for a loop_id, i.e. timestep -1
        bool chunk_learner = false;
        int init_method = 0;

        if (agent_type == 15) { // ChunkLearner: Init the data structure and calculate the chunk-sizes to try
#if (_RL_DEBUG > 1)
            std::cout << "[AgentProvider::search] Creating ChunkLearner for loop: " << loop_id << std::endl;
            std::cout << "[AgentProvider::search] Iterations: " << stats->n << ", Threads: " << stats->p << std::endl;
#endif
            chunk_sizes = new int[dimension];
            calculate_chunks(chunk_sizes, dimension, stats->n, stats->p);
            // Overwrites the agent type from 'ChunkLearner' to the new subtype to be used.
            read_env_int("KMP_RL_CHUNK_TYPE", agent_type);

#if (_RL_DEBUG > 1)
            std::cout << "[AgentProvider::search] ChunkLearner agent created. Chunk sizes: ";
            for (int i = 0; i < dimension; i++) {
                std::cout << chunk_sizes[i] << ", ";
            }
            std::cout << std::endl;
#endif

            if (agent_type == 15) {
                std::cout << "[AgentProvider::search] Aborting mission: Create nested Chunk-Learner" << std::endl;
                exit(-1);
            }

            chunk_learner = true;
        }
#if (_RL_DEBUG > 1)
        std::cout << "[AgentProvider::search] Creating agent for loop: " << loop_id << std::endl;
#endif
        // We insert -1 as value for the timestep, because the first time this gets checked is in the else block in timestep 0
        AgentProvider::get_timesteps().insert(std::make_pair(loop_id, 0));

        // The offset denotes the auto-methods already present in LB4OMP, so we can start our switch statement at 0
        auto *agent = create_agent(agent_type, loop_id, dimension, dimension);

        if (chunk_learner) {
            agent->set_name("Chunk-Learner (using " + agent->get_name() + ")");
        }

        AgentProvider::get_agents().insert(std::make_pair(loop_id, agent));
#if (_RL_DEBUG > 1)
        std::cout << "[AgentProvider::search] Agent created." << std::endl;
#endif
        print_agent_stats(agent, agent->get_current_state());

        if (chunk_learner) {
            print_agent_params(agent, chunk_sizes, dimension);
            return chunk_sizes[init_method];
        } else {
            print_agent_params(agent, nullptr, -1);
            return init_method;
        }
    } else { // Enters when agent exists for a loop_id, i.e. timestep >= 0
        AgentProvider::get_timesteps().at(loop_id) += 1;
//        std::cout << loop_id << ", ";
        int timestep = AgentProvider::get_timesteps().at(loop_id);
//        std::cout << " timestep: " << timestep << ", ";
        auto *agent = AgentProvider::get_agents().find(loop_id)->second;
//        std::cout << " name: " << agent->get_name() << ", ";
        int cMethod = agent->get_current_state();
//        std::cout << " cMethod: " << cMethod << ", ";
        int new_method = agent->step(0, timestep, stats);
//        std::cout << " new method: " << new_method;

// Update the best method if the current reward is better
        if (agent->get_current_reward() == agent->get_reward_num()[0]) {
            stats->bestTime = stats->cTime;
            stats->bestDLS = stats->cDLS;
            stats->bestLB = stats->cLB;
            stats->bestChunk = stats->cChunk;
            stats->bestStdDev = stats->cStdDev;
            stats->bestSkew = stats->cSkew;
            stats->bestKurt = stats->cKurt;
            stats->bestCOV = stats->cCOV;

#if (_RL_DEBUG > 0)
            std::cout << std::fixed << std::setprecision(6) << "[AgentProvider::search] Best reward updated. cTime: "
                      << stats->cTime << ", cDLS: " << stats->cDLS << ", cChunk: " << stats->cChunk << ", cLB: "
                      << stats->cLB << ", cCOV: " << stats->cCOV << std::endl;
#endif
        }

        print_agent_stats(agent, cMethod);

#if (_RL_DEBUG > 0)
        std::cout << "[AgentProvider::search] Timestep " << timestep << " completed. New method: " << new_method
                  << std::endl;
#endif

        if (agent_type == 15) {
            // ChunkLearner: Translate action index -> actual chunk-size and return it
#if (_RL_DEBUG > 1)
            std::cout << "[AgentProvider::search] Translating action index to chunk size: " << new_method << " --> "
                      << chunk_sizes[new_method] << std::endl;
#endif
            new_method = chunk_sizes[new_method];
        }

//        std::cout << "Current method: " << cMethod << ", New method: " << new_method << std::endl;
        return new_method;
    }
}

Agent *AgentProvider::create_agent(int agent_type, std::string loopname, int states, int actions) {
    Agent *agent = nullptr;
    int offset = 6;
    int new_type = agent_type - offset;

#if (_RL_DEBUG > 1)
    std::cout << "[AgentProvider::create_agent] New agent option: " << agent_type << " (offset: " << new_type << ")"
              << std::endl;
#endif
    switch (new_type) {
        case (0):
        case (1):
        case (2):
            agent = new QLearner(states, actions);
            break;
        case (3):
            agent = new DoubleQLearner(states, actions);
            break;
        case (4):
            agent = new SARSALearner(states, actions);
            break;
        case (5):
            agent = new ExpectedSARSALearner(states, actions);
            break;
        case (6):
            agent = new QVLearner(states, actions);
            break;
        case (7):
            agent = new RLearner(states, actions);
            break;
        case (8):
            agent = new DeepQLearner(states, actions);
            break;
        default:
            std::cout << "[AgentProvider::create_agent] Unknown agent type specified: " << agent_type
                      << ". Using default (Q-Learner)." << std::endl;
            agent = new QLearner(states, actions);
            break;
    }
    agent->set_loopname(loopname);
    agent->set_agent_type(agent_type);

    return agent;
}

std::unordered_map<std::string, int> &AgentProvider::get_timesteps() {
    return Get().timesteps;
}

std::unordered_map<std::string, Agent *> &AgentProvider::get_agents() {
    return Get().agents;
}

std::fstream &AgentProvider::get_filestream() {
    return Get().ofs;
}

int AgentProvider::calculate_chunks(int *array, int size, int n, int p) {
#if (_RL_DEBUG > 1)
    std::cout << "[AgentProvider::calculate_chunks] Calculating " << size << " chunks." << std::endl << "[";
#endif
    for (int i = 1; i <= size; i++) {
        int num = n / (pow(2, i) * p);
#if (_RL_DEBUG > 1)
        std::cout << num << ", ";
#endif
        array[i - 1] = num;
    }
#if (_RL_DEBUG > 1)
    std::cout << "]" << std::endl;
#endif
    return 0;
}

void AgentProvider::print_agent_params(Agent *agent, int *chunk_sizes_array, int dim) {
    std::ostringstream out;
    out << "[PARAMS-" << agent->get_loopname() << "]\n"
        << "Name = " << agent->get_name() << "\n"
        << "Rewards = " << agent->get_rewards_string() << "\n"
        << "RewardType = " << RewardLookup[agent->get_reward_type()] << "\n"
        << "InitType = " << InitLookup[agent->get_init_type()] << "\n"
        << "PolicyType = " << PolicyLookup[agent->get_policy_type()] << "\n"
        << "Alpha = " << std::fixed << std::setprecision(2) << agent->get_alpha_init() << "\n"
        << "AlphaMin = " << std::fixed << std::setprecision(2) << agent->get_alpha_min() << "\n"
        << "AlphaFactor = " << std::fixed << std::setprecision(2) << agent->get_alpha_decay() << "\n"
        << "Gamma = " << std::fixed << std::setprecision(2) << agent->get_gamma() << "\n"
        << "Epsilon = " << std::fixed << std::setprecision(2) << agent->get_epsilon_init() << "\n"
        << "EpsilonMin = " << std::fixed << std::setprecision(2) << agent->get_epsilon_min() << "\n"
        << "EpsilonFactor = " << std::fixed << std::setprecision(2) << agent->get_epsilon_decay() << "\n";

    std::string out_str = out.str();

    if (dim > 0) { // print the chunk-sizes to the ini file if using Chunk-Learner
        out_str += "Chunks = [";
        for (int i = 0; i < dim; i++) {
            out_str +=  std::to_string(chunk_sizes_array[i]);
            if (i != dim - 1) {
                out_str += ", ";
            }
        }
        out_str += "]\n";
    }

    out_str += "\n";
    std::cout << out_str;

    std::string fileData = agent->get_fileData();
    if (fileData != "notOutputtingStats") {
        std::fstream ofs;
        ofs.open(fileData + ".ini", std::fstream::in | std::fstream::out | std::fstream::app);
        ofs << out_str;
        ofs.close();
    }
}

void AgentProvider::print_agent_stats(Agent *agent, int cMethod) {
    std::string loopname = agent->get_loopname();

    int timestep = AgentProvider::get_timesteps().at(loopname);
    cMethod = agent->get_agent_type() == 15 ? chunk_sizes[cMethod] : cMethod;

    double alpha = agent->get_alpha();
    double epsilon = agent->get_epsilon();
    double reward = agent->get_current_reward();
    double low = agent->get_low();
    double high = agent->get_high();
    double reward_val = agent->get_current_reward_val();

    AgentProvider::get_filestream() << std::fixed << std::setprecision(6) << loopname << ",\t" << timestep << ",\t"
                                    << cMethod << ",\t" << alpha << ",\t" << epsilon << ",\t" << reward << ",\t" << low
                                    << ",\t" << high << ",\t" << reward_val << std::endl;
}