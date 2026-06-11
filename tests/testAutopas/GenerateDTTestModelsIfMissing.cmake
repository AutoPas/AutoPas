# Called by the runTests POST_BUILD step via cmake -P.
# Invokes the Python model-generation script only when the output files are absent,
# so repeated builds pay no cost beyond a file-existence check.
if (NOT EXISTS "${OUTPUT_DIR}/test_model.pkl" OR
    NOT EXISTS "${OUTPUT_DIR}/test_model_invalid.pkl")
    message(STATUS "Generating DecisionTreeTuning test models")
    execute_process(
        COMMAND "${PYTHON}" "${SCRIPT}" "${OUTPUT_DIR}"
        COMMAND_ERROR_IS_FATAL ANY
    )
endif ()