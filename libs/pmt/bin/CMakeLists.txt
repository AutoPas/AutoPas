project(PMT)

FetchContent_Declare(
  cxxopts
  GIT_REPOSITORY https://github.com/jarro2783/cxxopts.git
  GIT_TAG v3.2.1)
FetchContent_MakeAvailable(cxxopts)

add_executable(${PROJECT_NAME} main.cpp)

target_link_libraries(${PROJECT_NAME} cxxopts pmt)

target_include_directories(${PROJECT_NAME} PRIVATE ${CMAKE_BINARY_DIR})

install(TARGETS ${PROJECT_NAME} LIBRARY DESTINATION bin)
