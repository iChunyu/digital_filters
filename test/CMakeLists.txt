add_executable(test_biquad test_biquad.c)
target_link_libraries(test_biquad PRIVATE biquad_filter)
add_test(NAME test_biquad COMMAND test_biquad)

add_executable(test_sos test_sos.c)
target_link_libraries(test_sos PRIVATE biquad_filter)
add_test(NAME test_sos COMMAND test_sos)

add_executable(test_butter test_butter.c)
target_link_libraries(test_butter PRIVATE biquad_filter m)
add_test(NAME test_butter COMMAND test_butter)

add_executable(test_cheby test_cheby.c)
target_link_libraries(test_cheby PRIVATE biquad_filter m)
add_test(NAME test_cheby COMMAND test_cheby)
