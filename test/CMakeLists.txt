add_executable(test_biquad test_biquad.c)
target_link_libraries(test_biquad PRIVATE biquad_filter)
add_test(NAME test_biquad COMMAND test_biquad)

add_executable(test_sos test_sos.c)
target_link_libraries(test_sos PRIVATE biquad_filter)
add_test(NAME test_sos COMMAND test_sos)
