file(
        GLOB_RECURSE
        tuningLogToSQL_SRC
        "*.cpp"
        "*.h"
)

add_executable(tuningLogToSQL ${tuningLogToSQL_SRC})

target_include_directories(tuningLogToSQL PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})

target_link_libraries(tuningLogToSQL PUBLIC autopas autopasTools sqlite3)
