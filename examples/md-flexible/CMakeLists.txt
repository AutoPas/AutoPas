file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        "../../tests/testAutopas/testingHelpers/RandomGenerator.cpp"
        )

ADD_EXECUTABLE(md-flexible
        ${MY_SRC}
        )

TARGET_LINK_LIBRARIES(md-flexible
        autopas)

#-----------------test-----------------
# add check for current target
add_test(
        NAME md-flexible.test
        COMMAND md-flexible
                    --container linked
                    --cutoff 1.
                    --distribution-mean 5.0
                    --distribution-stddeviation 2.0
                    --data-layout soa
                    --functor lj
                    --iterations 10
                    --particles-generator gauss
                    --particles-per-dimension 10
                    --particle-spacing 0.4
                    --traversal c08,sliced
                    --verlet-rebuild-frequency 5
                    --verlet-skin-radius 0
        CONFIGURATIONS checkExamples
)

# add the executable to checkExamples as dependency
add_dependencies(checkExamples md-flexible)
