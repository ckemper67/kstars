add_executable(test_donuts test_donuts.cpp)
target_link_libraries(test_donuts donuts)
set_property(TARGET test_donuts PROPERTY CXX_STANDARD 17)
add_test(NAME DonutsLibTest COMMAND test_donuts)
