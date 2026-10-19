// -------------------------- Reinforcement Learning Extension ---------------------------------//
//  June 2022
//  Master Thesis
//  Luc Kury, <luc.kury@unibas.ch>
//  University of Basel, Switzerland
//  --------------------------------------------------------------------------------------------//

#include <unordered_map>

#include "kmp_explore_first_policy.h"
#include "../utils/utils.h"
#include "../agents/kmp_agent.h"

struct AgentVars {
    std::vector <std::vector<int>> visited;
    int noOnes = 0;
    int jumpIndex = 0;
    int lastTimestep = -1;
};

std::unordered_map <std::string, AgentVars> AgentVarsMap;

int ExploreFirstPolicy::policy(int episode, int timestep, Agent *agent) {
    int next_action;
    int no_actions = agent->get_action_space();
    int no_states = agent->get_state_space();
    int curr_state = agent->get_current_state();
    AgentVars &cData = AgentVarsMap[agent->get_loopname()];

    if (cData.visited.empty())
        cData.visited = std::vector < std::vector < int >> (no_actions, std::vector<int>(no_states, 0));

    if (timestep == cData.lastTimestep) // SARSA specific, as the next state is already updated
        curr_state = agent->get_next_state();

    if (cData.noOnes < no_states * no_actions) { // Try all actions possible in each state exactly once
#if (_RL_DEBUG > 0)
        if (cData.lastTimestep != timestep)
            std::cout << "[ExploreFirstPolicy::policy] Exploring! no_ones: " << cData.noOnes << " loop: " << agent->get_loopname() << std::endl;;
#endif

        for (next_action = 0; next_action < no_actions; next_action++) {
            if (cData.visited[curr_state][next_action] == 0) {
                int noOnesPerRow = 0; // Counter for the actions already visited in the next next state
                for (int j = 0; j < no_actions; j++) {
                    if (cData.visited[next_action][j] == 1)
                        noOnesPerRow++;
                }
                if (noOnesPerRow == no_actions) {
                    if (cData.lastTimestep != timestep) {
                        cData.jumpIndex += 1;
                    }
                    next_action = (next_action + cData.jumpIndex) % no_actions;
                }
                if (cData.lastTimestep != timestep) {
                    cData.visited[curr_state][next_action] = 1;
                    cData.noOnes += 1;
                }
                break;
            }
        }
    } else { // Determine which action have the best reward; use the averaged reward of this action in all states
#if (_RL_DEBUG > 0)
        if (cData.lastTimestep != timestep)
            std::cout << "[ExploreFirstPolicy::policy] Exploiting!" << std::endl;;
#endif
        next_action = argmax(*agent->get_table(), curr_state, no_actions);
    }
    cData.lastTimestep = timestep;
    return next_action;
}
