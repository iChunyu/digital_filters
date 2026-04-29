add_executable(test_biquad test_biquad.c)
target_link_libraries(test_biquad PRIVATE biquad_filter)
add_test(NAME test_biquad COMMAND test_biquad)
