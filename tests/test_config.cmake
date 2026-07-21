include(CTest)

find_program(PYTHON_FROM_PATH python REQUIRED)
set(Python3_EXECUTABLE "${PYTHON_FROM_PATH}" CACHE FILEPATH "Python interpreter" FORCE)
find_package(Python3 REQUIRED COMPONENTS Interpreter)
message(STATUS "Python3 executable: ${Python3_EXECUTABLE}")

file(GLOB TEST_CONFIGS CONFIGURE_DEPENDS
    "${CMAKE_CURRENT_LIST_DIR}/*.config"
)

foreach(CONFIG_FILE ${TEST_CONFIGS})
    get_filename_component(TEST_NAME ${CONFIG_FILE} NAME_WE)
    add_test(
        NAME "hubbard_${TEST_NAME}"
        COMMAND
            ${Python3_EXECUTABLE}
            "${CMAKE_CURRENT_LIST_DIR}/run_test.py"
            --exe "$<TARGET_FILE:Hubbard>"
            --config "${CONFIG_FILE}"
            --plot-script "${CMAKE_CURRENT_LIST_DIR}/plot.py"
    )

    set_tests_properties("hubbard_${TEST_NAME}"
        PROPERTIES
            WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}"
    )
endforeach()



execute_process(
    COMMAND "${Python3_EXECUTABLE}" -c "import matplotlib; import numpy; import mrock"
    RESULT_VARIABLE PYTHON_MODULE_CHECK
    OUTPUT_VARIABLE PYTHON_MODULE_CHECK_OUT
    ERROR_VARIABLE PYTHON_MODULE_CHECK_ERR
)

if(NOT PYTHON_MODULE_CHECK EQUAL 0)
    message(FATAL_ERROR
        "The selected Python interpreter does not have the required modules.\n"
        "Python: ${Python3_EXECUTABLE}\n"
        "Error:\n${PYTHON_MODULE_CHECK_ERR}"
    )
endif()