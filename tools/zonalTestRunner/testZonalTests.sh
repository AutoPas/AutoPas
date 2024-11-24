
# Build and run all relevant tests for this thesis
# Run this from the build folder

cmake .. -DCMAKE_EXPORT_COMPILE_COMMANDS=ON -DMD_FLEXIBLE_USE_MPI=ON
make mdFlexTests -j 5

ctest -R mdFlexTests/ZonalMethodTest
ctest -R mdFlexTests/HalfShellTest
ctest -R mdFlexTests/FullShellTest
ctest -R mdFlexTests/RegionTest
ctest -R mdFlexTests/TestHaloZonalParticles/RegularGridDecompositionTest.
