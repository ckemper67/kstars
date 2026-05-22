add_library(donuts STATIC donuts.cpp)
target_include_directories(donuts PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
set_property(TARGET donuts PROPERTY CXX_STANDARD 17)

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
