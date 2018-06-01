file(GLOB_RECURSE MY_SRC
        "*.cpp"
        "*.h"
        )

ADD_EXECUTABLE(md-flexible
        ${MY_SRC}
        )

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -g")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -S")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info-vec-all")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info-vec-missed")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info-vec-optimized")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fopt-info-vec-note")

#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-report=2")
#set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -qopt-report-phase=vec")
TARGET_LINK_LIBRARIES(md-flexible
        autopas)

#-----------------test-----------------
# add check for current target
add_test(
        NAME md-flexible.test
        COMMAND md-flexible --container linked --cutoff 1. --data-layout soa --functor lj --iterations 10
            --particles-per-dimension 10 --particle-spacing 0.4
        CONFIGURATIONS checkExamples
)

# add the executable to checkExamples as dependency
add_dependencies(checkExamples md-flexible)
