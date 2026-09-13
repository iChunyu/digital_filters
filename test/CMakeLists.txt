add_executable(test_biquad test_biquad.c)
target_link_libraries(test_biquad PRIVATE digital_filters)
add_test(NAME test_biquad COMMAND test_biquad)

add_executable(test_butter test_butter.c)
target_link_libraries(test_butter PRIVATE digital_filters m)
add_test(NAME test_butter COMMAND test_butter)

add_executable(test_cheby test_cheby.c)
target_link_libraries(test_cheby PRIVATE digital_filters m)
add_test(NAME test_cheby COMMAND test_cheby)

add_executable(test_notch test_notch.c)
target_link_libraries(test_notch PRIVATE digital_filters m)
add_test(NAME test_notch COMMAND test_notch)

# 与 Python 对比用的可执行文件：生成器把 CSV 写进 build 目录；
# scipy 对比从**同一个**目录读回，这样源码目录里的陈旧 CSV
# 永远不会遮蔽新跑出来的结果。
add_executable(test_butter_with_py test_butter_with_py.c)
target_link_libraries(test_butter_with_py PRIVATE digital_filters m)

add_executable(test_cheby_with_py test_cheby_with_py.c)
target_link_libraries(test_cheby_with_py PRIVATE digital_filters m)

add_executable(test_notch_with_py test_notch_with_py.c)
target_link_libraries(test_notch_with_py PRIVATE digital_filters m)

# 把配套的 Python 脚本拷到可执行文件旁边
configure_file(test_butter_use_py.py ${CMAKE_CURRENT_BINARY_DIR}/test_butter_use_py.py COPYONLY)
configure_file(test_cheby_use_py.py ${CMAKE_CURRENT_BINARY_DIR}/test_cheby_use_py.py COPYONLY)
configure_file(test_notch_use_py.py ${CMAKE_CURRENT_BINARY_DIR}/test_notch_use_py.py COPYONLY)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
    add_test(NAME butter_csv_gen
             COMMAND test_butter_with_py ${CMAKE_CURRENT_BINARY_DIR}/test_butter_data.csv)
    add_test(NAME cheby_csv_gen
             COMMAND test_cheby_with_py ${CMAKE_CURRENT_BINARY_DIR}/test_cheby_data.csv)
    add_test(NAME notch_csv_gen
             COMMAND test_notch_with_py ${CMAKE_CURRENT_BINARY_DIR}/test_notch_data.csv)
    add_test(NAME scipy_compare
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/compare_scipy.py
                     ${CMAKE_CURRENT_BINARY_DIR})
    # numpy/scipy 不可用时 compare_scipy.py 以 77 退出
    set_tests_properties(scipy_compare PROPERTIES
        DEPENDS "butter_csv_gen;cheby_csv_gen"
        SKIP_RETURN_CODE 77)

    # 陷波器的频响对拍：scipy.signal.bilinear 是 c2d 的独立实现
    #（已确认与本库逐位吻合），C 侧则是 f32 时域递推。
    add_test(NAME notch_py_compare
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_notch_use_py.py
                     ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(notch_py_compare PROPERTIES
        DEPENDS "notch_csv_gen"
        SKIP_RETURN_CODE 77)

    # 独立复现 zpk 增益管线的 f64-vs-f32 精度对比。只需要 numpy
    #（不需要 scipy）；缺失时以 77 跳过。
    add_test(NAME verify_zpk_gain
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/verify_zpk_gain.py)
    set_tests_properties(verify_zpk_gain PROPERTIES SKIP_RETURN_CODE 77)
endif()
