export SIMULATION=FallingDrop
export STRATEGY=ReinforcementLearning

for i in 1 2 3 4 5 6 7 8 9 10 11; do
  (
    sbatch "run_alpha_gamma_cluster_${i}.cmd"
  )
done
