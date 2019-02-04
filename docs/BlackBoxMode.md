## BlackBoxMode
When the blackbox mode is enabled, the particle exchange is independent
of the container type. It is required to be as follows:
1. Leaving particles should be sent out and be received in every
iteration.
2. Halo particles need to be copied in each iteration.
3. In each iteration the halo particles should be deleted after the
force is calculated.

The blackbox mode enables fluent interplay across different container
types.
