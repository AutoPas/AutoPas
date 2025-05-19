include(FetchContent REQUIRED)

if(NOT DEFINED NLOHMANN_JSON_GIT_TAG)
    set(NLOHMANN_JSON_GIT_TAG "master")
endif()
message(STATUS "FetchContent nlohmann_json")
FetchContent_Declare(
        nlohmann_json
        GIT_REPOSITORY https://github.com/nlohmann/json.git
        GIT_TAG        ${NLOHMANN_JSON_GIT_TAG}
)
FetchContent_MakeAvailable(nlohmann_json)