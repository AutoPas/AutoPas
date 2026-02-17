# AutoTuning

One of the primary features of AutoPas is the dynamic automated algorithm selection for the force calculation.
Since AutoPas follows the philosophy of [Convention over Configuration](https://en.wikipedia.org/wiki/Convention_over_configuration), it will automatically tune itself to find the optimum of all its internal algorithm configurations.
However, this naive tuning is just an exhaustive search by trialling all algorithms over successive iterations.
While this is guaranteed to find the optimum, it is very costly and inefficient.
Therefore, AutoPas offers a bunch of [Tuning Strategies](https://github.com/AutoPas/AutoPas/blob/master/src/autopas/options/TuningStrategyOption.h) that can be applied and combined for an optimal experience.

## Tuning Loop
To better understand what the strategies do, we first take a look at the way AutoPas does its tuning.

The basic idea is that the [AutoTuner](https://github.com/AutoPas/AutoPas/blob/master/src/autopas/tuning/AutoTuner.h) class determines a configuration, which the [LogicHandler](https://github.com/AutoPas/AutoPas/blob/master/src/autopas/LogicHandler.h) then assembles and applies to the [ParticleContainer](https://github.com/AutoPas/AutoPas/blob/master/src/autopas/containers/ParticleContainerInterface.h).
In order to determine the optimal configuration, the AutoTuner has a queue of candidates that it tries out by using it for one interaction computation in order to measure their performance.
Those iterations while candidate configurations are evaluated are referred to as a tuning phase.
AutoPas works on the assumption that particle simulations evolve very slow and thus the simulation state of subsequent iterations is sufficiently similar to compare the configurations' performances.
Therefore, the simulation can advance with each iteration during tuning and doesn't have to pause.
When the queue is fully processed, the configuration with the optimal measurement is determined to be the current algorithm optimum and is used until the start of the next tuning phase.

### Tuning Interval
The interval between the start of two tuning phases can be configured by the user.
Iterations for sampling the configurations do count against this interval, so if a tuning interval of 1000 is set, tuning will start at iterations 0, 1000, and 2000.
Should the tuning take longer than one interval, the subsequent tuning phases are skipped until the first is complete.

### Tuning Strategies
By default, the AutoTuner configuration queue is filled with all applicable configurations and processed sequentially.
Tuning strategies can filter, reorder, or otherwise modify this queue.
After each sampling of a configuration, tuning strategies receive the measurement, called evidence, which they can use to change how they modify the queue. 
For example, the predictive tuning strategy can use evidence from past tuning phases to filter the current tuning phase's queue.
Tuning strategies can also use live information about the current state of the domain, such as density or homogeneity, to filter the queue. 
For example, the rules-based tuning strategy can be used to filter out the Direct Sum container if the simulation has enough particles.
Multiple strategies can be activated at the same time.
They are all applied in every decision process in the order they are specified in the vector of `TuningStrategyOption`s passed to AutoPas.

### Sampling, Selecting, and Smoothing
For each evidence, AutoPas combines multiple subsequent measurements, called samples.
This is necessary for two reasons.
First, Verlet based algorithms require iterations with rebuild steps that take significantly longer than those without.
Thus at least two iterations are necessary to calculate a weighted average.
Second, especially measurements from iterations, which take significantly less time than one second, tend to be noisy.
This is tackled by reducing all measurements of a type (rebuilding vs. not-rebuilding) for one configuration with a [SelectorStrategy](https://github.com/AutoPas/AutoPas/blob/master/src/autopas/options/SelectorStrategyOption.h).
Whilst gathering many samples can help reduce the noise, it can also be costly due to having to gather many samples for sub-optimal configurations.
To alleviate this, when data from at least two previous tuning phases already exists, a variation of [LOESS](https://en.wikipedia.org/wiki/Local_regression) can be applied to mitigate positive outliers. 
Use this feature with caution since it can sometimes lead to undesirable tuning outcomes. 
When in doubt, compare the reduced vs smoothed values from the `TuningDataLogger`. 

## Tuning Metric
AutoPas not only offers tuning for minimal runtime.
The tuning metric can also be switched to minimize energy consumption.
See [Building](https://github.com/AutoPas/AutoPas/blob/master/docs/userdoc/Building.md) for instructions on how to compile support for this.

## Multiple Interaction Types
Since PR [751](https://github.com/AutoPas/AutoPas/pull/751), AutoPas also supports algorithms for 3-body interactions.
Oftentimes, these are used additionally to standard pairwise interactions. For such simulations, AutoPas runs two separate AutoTuners, one for pairwise and one for triwise interaction.
Both AutoTuners will start their tuning phases at the same iterations by setting the start of the tuning interval when both AutoTuners finish tuning over their search space.  
_Note: Performance may vary if the AutoTuners end up with configurations that use a different [container option](https://github.com/AutoPas/AutoPas/blob/tuningDoc/src/autopas/options/ContainerOption.h). 
In this case, AutoPas has to rebuild neighbor lists every time the interaction type changes._
