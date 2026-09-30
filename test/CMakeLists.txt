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

# test_peak.c 另外包含 notch_filter.h——对偶不变量测试要拿同一组参数的
# 陷波系数来对照（H_peak(ξ,g) ≡ 1/H_notch(ξ,1/g)）。
add_executable(test_peak test_peak.c)
target_link_libraries(test_peak PRIVATE digital_filters m)
add_test(NAME test_peak COMMAND test_peak)

# 与 Python 对比用的可执行文件：生成器把 CSV 写进 build 目录；
# scipy 对比从**同一个**目录读回，这样源码目录里的陈旧 CSV
# 永远不会遮蔽新跑出来的结果。
add_executable(test_butter_with_py test_butter_with_py.c)
target_link_libraries(test_butter_with_py PRIVATE digital_filters m)

add_executable(test_cheby_with_py test_cheby_with_py.c)
target_link_libraries(test_cheby_with_py PRIVATE digital_filters m)

add_executable(test_notch_with_py test_notch_with_py.c)
target_link_libraries(test_notch_with_py PRIVATE digital_filters m)

add_executable(test_peak_with_py test_peak_with_py.c)
target_link_libraries(test_peak_with_py PRIVATE digital_filters m)

# 把配套的 Python 脚本拷到可执行文件旁边
configure_file(test_butter_use_py.py ${CMAKE_CURRENT_BINARY_DIR}/test_butter_use_py.py COPYONLY)
configure_file(test_cheby_use_py.py ${CMAKE_CURRENT_BINARY_DIR}/test_cheby_use_py.py COPYONLY)
configure_file(test_notch_use_py.py ${CMAKE_CURRENT_BINARY_DIR}/test_notch_use_py.py COPYONLY)
configure_file(test_peak_use_py.py ${CMAKE_CURRENT_BINARY_DIR}/test_peak_use_py.py COPYONLY)

find_package(Python3 COMPONENTS Interpreter)
if(Python3_FOUND)
    add_test(NAME butter_csv_gen
             COMMAND test_butter_with_py ${CMAKE_CURRENT_BINARY_DIR}/test_butter_data.csv)
    add_test(NAME cheby_csv_gen
             COMMAND test_cheby_with_py ${CMAKE_CURRENT_BINARY_DIR}/test_cheby_data.csv)
    add_test(NAME notch_csv_gen
             COMMAND test_notch_with_py ${CMAKE_CURRENT_BINARY_DIR}/test_notch_data.csv)
    add_test(NAME peak_csv_gen
             COMMAND test_peak_with_py ${CMAKE_CURRENT_BINARY_DIR}/test_peak_data.csv)
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

    # 峰值滤波器的频响对拍，与上面同构。被测原型与陷波互为倒数
    # （H1(ξ,g) ≡ 1/H2(ξ,1/g)），所以两个 CSV 也可以直接逐点互倒比较。
    add_test(NAME peak_py_compare
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/test_peak_use_py.py
                     ${CMAKE_CURRENT_BINARY_DIR})
    set_tests_properties(peak_py_compare PROPERTIES
        DEPENDS "peak_csv_gen"
        SKIP_RETURN_CODE 77)

    # 独立复现 zpk 增益管线的 f64-vs-f32 精度对比。只需要 numpy
    #（不需要 scipy）；缺失时以 77 跳过。
    add_test(NAME verify_zpk_gain
             COMMAND ${Python3_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/verify_zpk_gain.py)
    set_tests_properties(verify_zpk_gain PROPERTIES SKIP_RETURN_CODE 77)
else()
    # 静默少注册 8 项测试会让“13/13”变成“5/5”而看不出来。这里把话说清楚，
    # 不让“未运行”与“全通过”在构建输出里无从区分。
    message(WARNING
        "Python3 未找到：8 项测试不会注册（4 个 CSV 生成器 + scipy_compare + "
        "notch_py_compare + peak_py_compare + verify_zpk_gain）。"
        "本次 ctest 只会跑 5 项 C 测试——不要把它当成 13/13。")
endif()
