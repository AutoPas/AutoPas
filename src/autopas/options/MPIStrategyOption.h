/**
 * @file MPIStrategyOption.h
 * @author W.Thieme
 * @date 28.05.2020
*/

#pragma once

#include "Option.h"

namespace autopas {
inline namespace options {

class MPIStrategyOption : public Option<MPIStrategyOption> {
public:
  /**
   * Possible options for the use of mpi in tuning
   */
  enum Value {
    /**
     * Do not use mpi. Every AutoPas instance acts independently
     */
            noMPI,
    /**
     * local tuning in a reduced search space in each rank with subsequent global comparison for the best strategy
     */
            divideAndConquer,
  };

  /**
   * Constructor
   */
  MPIStrategyOption() = default;

  /**
   * Constructor from value.
   * @param option
   */
  constexpr MPIStrategyOption(Value option) : _value(option) {}

  /**
   * Cast to value.
   * @return
   */
  constexpr operator Value() const { return _value; }

  /**
   * Provides a way to iterate over the possible choices of TuningStrategy.
   * @return map option -> string representation
   */
  static std::map<MPIStrategyOption, std::string> getOptionNames() {
    return {
            {MPIStrategyOption::noMPI,            "no-mpi"},
            {MPIStrategyOption::divideAndConquer, "divide-and-conquer"},
    };
  }

private:
  Value _value{Value(-1)};

};
} // namespace options
} // namespace autopas