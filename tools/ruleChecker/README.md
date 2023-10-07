# Rule-based Tuning Rule Checker
The rule checker for rule-based tuning can be used to verify that the partial order of configurations defined in a rule file does not contradict the actual runtime behavior of a simulation. Here is how it fits into the larger process of writing good rule files:

1. Run a suite of benchmark scenarios using full-search tuning with enabled data collection (`md-flexible --use-tuning-log true`) and one tuning phase.  
Each simulation gives you a `tuningLog.txt` file.
2. Come up with new rules (through theoretical modeling, developer knowledge, benchmarks etc.) in a file `tuningRules.txt`.
3. Call `ruleChecker tuningRules.rule tuningLog.txt tuningLog2.txt ...`  
This will replay all logs from the benchmarks and use the rule based tuning strategy in verify mode on them. This means that it executes the rule file, creates the partial order of configurations, but instead of filtering out configurations, it checks in the logs whether the partial order contradicts the observation in the benchmark.  
For example, the rule file might state that configurations with the `VerletClusterLists` container are always faster than `LinkedCells` containers in scenarios with less than 2 particles per cell, but in the benchmarks, `VerletClusterLists` using the `AoS` data layout and no Newton3 was actually slightly slower than `LinkedCells`, `c08`, `SoA`, `Newton3`.  
The rule checker logs all such observations together with the actual difference in performance. Now, you can use this information to adjust the rule file.
4. Adjust the rules using the information about wrong partial orders.
5. If happy, stop. Otherwise, go to 3.

In addition to any wrong partial orders, the rule checker also prints how many configurations would have been filtered out in the simulation by this rule file, and how much tuning time would have been saved. This allows you to judge whether your rules file is strict enough to give meaningful performance benefits.

