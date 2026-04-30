add_executable(test_biquad test_biquad.c)
target_link_libraries(test_biquad PRIVATE digital_filters)
add_test(NAME test_biquad COMMAND test_biquad)

add_executable(test_sos test_sos.c)
target_link_libraries(test_sos PRIVATE digital_filters)
add_test(NAME test_sos COMMAND test_sos)

add_executable(test_butter test_butter.c)
target_link_libraries(test_butter PRIVATE digital_filters m)
add_test(NAME test_butter COMMAND test_butter)

add_executable(test_cheby test_cheby.c)
target_link_libraries(test_cheby PRIVATE digital_filters m)
add_test(NAME test_cheby COMMAND test_cheby)

add_executable(test_static_butter test_static_butter.c)
target_link_libraries(test_static_butter PRIVATE digital_filters m)
add_test(NAME test_static_butter COMMAND test_static_butter)

add_executable(test_static_cheby test_static_cheby.c)
target_link_libraries(test_static_cheby PRIVATE digital_filters m)
add_test(NAME test_static_cheby COMMAND test_static_cheby)
