cmake_minimum_required(VERSION 3.4.1)
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
project(proj_permute VERSION 1.0
               DESCRIPTION "Wavelet Trees with Learned Indexes for Sortedness"
               LANGUAGES CXX)

add_compile_options(-march=native)

set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
set(CMAKE_MODULE_PATH "${PROJECT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH})

option(DEBUG "Debug mode on or off" OFF)
option(PROFILE "Profile mode on or off" OFF)
if(DEBUG)
    set(CMAKE_BUILD_TYPE "Debug")
    set(CMAKE_CXX_FLAGS "-g -O0")
elseif(PROFILE)
    set(CMAKE_BUILD_TYPE "Debug")
    set(CMAKE_CXX_FLAGS "-g -O2 -DPROFILE")
else()
    set(CMAKE_BUILD_TYPE "Release")
    set(CMAKE_CXX_FLAGS "-O3 -DVALIDATE")
endif()

find_package(Git QUIET)
if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
# Update submodules as needed
    option(GIT_SUBMODULE "Check submodules during build" ON)
    if(GIT_SUBMODULE)
        message(STATUS "Submodule update")
        execute_process(COMMAND ${GIT_EXECUTABLE} submodule update --init --recursive
                        WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
                        RESULT_VARIABLE GIT_SUBMOD_RESULT)
        if(NOT GIT_SUBMOD_RESULT EQUAL "0")
            message(FATAL_ERROR "git submodule update --init --recursive failed with ${GIT_SUBMOD_RESULT}, please checkout submodules")
        endif()
    endif()
endif()

if(NOT EXISTS "${PROJECT_SOURCE_DIR}/external/CMakeLists.txt")
    message(FATAL_ERROR "Utilities submodules were not downloaded! GIT_SUBMODULE was turned off or failed. Please update submodules and try again.")
endif()

add_subdirectory(external)



file(GLOB_RECURSE
        TMPDB_SRCS
        CONFIGURE_DEPENDS
        ${CMAKE_SOURCE_DIR}/src/*.cpp
        ${CMAKE_SOURCE_DIR}/src/*.hpp
)

add_library(permute OBJECT "${TMPDB_SRCS}")

target_compile_features(permute PUBLIC
    cxx_std_20
)

target_link_libraries(permute PUBLIC
    spdlog::spdlog
    roaring
    radix_spline_lib
    learned-secondary-index
)

target_compile_options(permute PUBLIC 
$<$<BOOL:${DEBUG}>:-g>)


target_include_directories(permute PUBLIC
    ${CMAKE_SOURCE_DIR}/src
)

add_executable(permute_bench compare_binsearch_wt.cpp
)
target_compile_definitions(permute_bench PUBLIC VALIDATE)

target_compile_options(permute_bench PUBLIC 
$<$<BOOL:${DEBUG}>:-g>)


target_link_libraries(permute_bench PUBLIC
    permute
)

add_executable(permute_test_bvrs permute_test.cpp)
target_compile_definitions(permute_test_bvrs PUBLIC BVRS)
target_compile_options(permute_test_bvrs PUBLIC 
$<$<BOOL:${DEBUG}>:-g>
$<$<BOOL:${PROFILE}>:-g -O2 -DPROFILE>)
target_link_libraries(permute_test_bvrs PUBLIC
    permute
)

add_executable(permute_test_roaring permute_test.cpp)
target_compile_definitions(permute_test_roaring PUBLIC ROARING)
target_compile_options(permute_test_roaring PUBLIC 
$<$<BOOL:${DEBUG}>:-g>
$<$<BOOL:${PROFILE}>:-g -O2 -DPROFILE>)
target_link_libraries(permute_test_roaring PUBLIC
    permute
)

add_executable(permute_flat_multiary permute_test.cpp)
target_compile_definitions(permute_flat_multiary PUBLIC FLAT_MULTIARY)
target_compile_options(permute_flat_multiary PUBLIC 
$<$<BOOL:${DEBUG}>:-g>
$<$<BOOL:${PROFILE}>:-g -O3 -DPROFILE>)
target_link_libraries(permute_flat_multiary PUBLIC
    permute
)

add_executable(li_wavelet_test li_wavelet_test.cpp)
target_compile_definitions(li_wavelet_test PUBLIC
$<$<BOOL:${DEBUG}>:-g>
$<$<BOOL:${PROFILE}>:-g -O3 -DPROFILE>
)
target_link_libraries(li_wavelet_test PUBLIC
    permute
    learned-secondary-index
    learned-hashing
)

add_executable(li_lsi_test li_lsi_test.cpp)
target_compile_options(li_lsi_test PRIVATE -march=native -O3)
target_compile_definitions(li_lsi_test PUBLIC
    $<$<BOOL:${DEBUG}>:-g>
    $<$<BOOL:${PROFILE}>:-DPROFILE>
)

target_link_libraries(li_lsi_test PUBLIC
    learned-secondary-index
    learned-hashing
    spdlog::spdlog
)
