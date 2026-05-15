add_executable(test_donuts test_donuts.cpp)
target_link_libraries(test_donuts PRIVATE donuts)
target_compile_features(test_donuts PRIVATE cxx_std_17)

add_test(NAME DonutsAlgorithm COMMAND test_donuts)
