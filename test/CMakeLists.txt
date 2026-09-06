add_executable(test_biquad test_biquad.c)
target_link_libraries(test_biquad PRIVATE digital_filters)
add_test(NAME test_biquad COMMAND test_biquad)

add_executable(test_butter test_butter.c)
target_link_libraries(test_butter PRIVATE digital_filters m)
add_test(NAME test_butter COMMAND test_butter)

add_executable(test_cheby test_cheby.c)
target_link_libraries(test_cheby PRIVATE digital_filters m)
add_test(NAME test_cheby COMMAND test_cheby)

# Comparison-with-Python executables: generators write CSVs into the build
# dir; the scipy comparison reads them back from the SAME dir so a stale
# source-dir CSV can never shadow a fresh run.
add_executable(test_butter_with_py test_butter_with_py.c)
target_link_libraries(test_butter_with_py PRIVATE digital_filters m)

add_executable(test_cheby_with_py test_cheby_with_py.c)
target_link_libraries(test_cheby_with_py PRIVATE digital_filters m)

# Copy Python companion scripts alongside the executables
configure_file(test_butter_use_py.py ${CMAKE_CURRENT_BINARY_DIR}/test_butter_use_py.py COPYONLY)
configure_file(test_cheby_use_py.py ${CMAKE_CURRENT_BINARY_DIR}/test_cheby_use_py.py COPYONLY)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
    add_test(NAME butter_csv_gen
             COMMAND test_butter_with_py ${CMAKE_CURRENT_BINARY_DIR}/test_butter_data.csv)
    add_test(NAME cheby_csv_gen
             COMMAND test_cheby_with_py ${CMAKE_CURRENT_BINARY_DIR}/test_cheby_data.csv)
    add_test(NAME scipy_compare
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare_scipy.py
                     ${CMAKE_CURRENT_BINARY_DIR})
    # compare_scipy.py exits 77 when numpy/scipy is unavailable
    set_tests_properties(scipy_compare PROPERTIES
        DEPENDS "butter_csv_gen;cheby_csv_gen"
        SKIP_RETURN_CODE 77)
endif()
