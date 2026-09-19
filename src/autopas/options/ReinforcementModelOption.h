/**
 * @file ReinforcementModelOption.h
 * @author P. Metscher
 * @date 02.01.2026
 */

#pragma once

#include <set>

#include "autopas/options/Option.h"

namespace autopas {
inline namespace options {
/**
 * Option for the reinforcement learning model used in the Deep Reinforcement Learning tuning strategy.
 */
class ReinforcementModelOption : public Option<ReinforcementModelOption> {
 public:
  /**
   * Possible reinforcement learning models for the Reinforcement Learning tuning strategy.
   */
  enum Value {
    /**
     * The SARSA model.
     *
     * This model updates the state using the SARSA update rule
     *
     * @f[Q(s, a) \leftarrow Q(s, a) + \alpha \cdot (r + \gamma \cdot Q(s', a') - Q(s, a))@f]
     *
     * where @f$Q(s, a)@f$ is the action-value function, @f$\alpha@f$ is the learning rate, @f$r@f$ is the reward,
     * @f$\gamma@f$ is the discount factor, @f$s'@f$ is the next state and @f$a'@f$ is the next action.
     */
    SARSA,

    /**
     * The Q-learning model.
     *
     * This model updates the state using the Q-learning update rule
     *
     * @f[Q(s, a) \leftarrow Q(s, a) + \alpha \cdot (r + \gamma \cdot \max_{a'} Q(s', a') - Q(s, a))@f]
     *
     * where @f$Q(s, a)@f$ is the action-value function, @f$\alpha@f$ is the learning rate, @f$r@f$ is the reward,
     * @f$\gamma@f$ is the discount factor, @f$s'@f$ is the next state and @f$a'@f$ is the next action.
     */
    QLearning,
  };

  /**
   * Constructor.
   */
  ReinforcementModelOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr ReinforcementModelOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Provides a way to iterate over the possible choices of ReinforcementModels.
   * @return map option -> string representation
   */
  static std::map<ReinforcementModelOption, std::string> getOptionNames() {
    return {
        {Value::SARSA, "SARSA"},
        {Value::QLearning, "QLearning"},
    };
  };

 private:
  Value _value{Value(-1)};
};
}  // namespace options
}  // namespace autopas
