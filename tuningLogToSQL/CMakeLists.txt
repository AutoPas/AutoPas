find_package(SQLite3)
if(SQLite3_FOUND)
    message("Adding target 'tuningLogToSQL'")
    file(
            GLOB_RECURSE
            tuningLogToSQL_SRC
            "*.cpp"
            "*.h"
    )
    add_executable(tuningLogToSQL ${tuningLogToSQL_SRC})
    target_include_directories(tuningLogToSQL PUBLIC ${CMAKE_CURRENT_SOURCE_DIR} ${SQLite3_INCLUDE_DIRS})
    target_link_libraries(tuningLogToSQL PUBLIC autopas autopasTools ${SQLite3_LIBRARIES})

else()
    message("SQLite3 dev package not found. Try installing libsqlite3-dev.")
    message("Not adding 'tuningLogToSQL' as a target")
endif()
