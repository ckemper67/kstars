add_library(donuts STATIC donuts.cpp simd_warp.cpp)

target_include_directories(donuts
    PUBLIC  ${CMAKE_CURRENT_SOURCE_DIR}
    PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/../internalguide   # pocketfft_hdronly.h
)

target_compile_features(donuts PUBLIC cxx_std_17)

# Standalone build: user can override pocketfft location with -DPOCKETFFT_DIR=...
if(DEFINED POCKETFFT_DIR)
    target_include_directories(donuts PRIVATE ${POCKETFFT_DIR})
endif()

if(BUILD_TESTING)
    add_subdirectory(tests)
endif()
