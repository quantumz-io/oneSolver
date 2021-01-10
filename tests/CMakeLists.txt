find_package(Boost 1.65 REQUIRED COMPONENTS unit_test_framework)

enable_testing()

set(TEST_APP_NAME "main-test")

set(TEST_MODULE_PATH "${PROJECT_SOURCE_DIR}/tests")
set(TEST_SOURCE_PATH "${TEST_MODULE_PATH}")

set (TEST_SOURCE_FILES "${TEST_SOURCE_PATH}/qubo_test.cpp"
                       "${TEST_SOURCE_PATH}/qubo_helpers_test.cpp"                     
                       "${TEST_SOURCE_PATH}/io_test.cpp"
                       "${TEST_SOURCE_PATH}/devices_test.cpp"
                       "${TEST_SOURCE_PATH}/exhaustive_test.cpp"
                       )

add_executable(${TEST_APP_NAME} ${TEST_SOURCE_FILES})

target_include_directories(${TEST_APP_NAME} PRIVATE "${PROJECT_SOURCE_DIR}/include"
                                                    ${Boost_INCLUDE_DIRS})

target_compile_definitions(${TEST_APP_NAME} PRIVATE "BOOST_TEST_DYN_LINK=1")

target_link_libraries(${TEST_APP_NAME} ${Boost_UNIT_TEST_FRAMEWORK_LIBRARY})

add_test(NAME qubo_test COMMAND ${TEST_APP_NAME} --run_test=qubo_test)
add_test(NAME qubo_read_from_file COMMAND ${TEST_APP_NAME} --run_test=read_file_test)
add_test(NAME qubo_helpers_test COMMAND ${TEST_APP_NAME} --run_test=qubo_helpers_test)
add_test(NAME io_test COMMAND ${TEST_APP_NAME} --run_test=io_test)
add_test(NAME exhaustive_test COMMAND ${TEST_APP_NAME} --run_test=exhaustive_test)
add_test(NAME devices_test COMMAND ${TEST_APP_NAME} --run_test=devices_test)
