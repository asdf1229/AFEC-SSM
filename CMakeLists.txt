cmake_minimum_required(VERSION 3.10)

project(SSM_GED LANGUAGES CXX)

# --------------------------------------------
# Optimization flags
# --------------------------------------------
# add_compile_options(-O3 -g -Wall -march=native)
add_compile_options(-O3 -g -Wall)
# if (CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang")
#    add_compile_options(-O3 -g -Wall -march=native -fsanitize=address -fno-omit-frame-pointer)
#    set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -fsanitize=address")
# endif()

include_directories(
    ${PROJECT_SOURCE_DIR}
)

# --------------------------------------------
# File sources
# --------------------------------------------
file(GLOB_RECURSE SSM_GED_SOURCES
    ${PROJECT_SOURCE_DIR}/graph/*.cpp
    ${PROJECT_SOURCE_DIR}/matching/*.cpp
)

find_package(Threads REQUIRED)

# --------------------------------------------
# Macro option for switching matching algorithm
# --------------------------------------------
option(APPROXIMATE_MATCHING_V2 "Enable matching V2 implementation" OFF)

if (APPROXIMATE_MATCHING_V2)
    message(STATUS "Using Approximate_Matching_V2")
    add_definitions(-DAPPROXIMATE_MATCHING_V2)
else()
    message(STATUS "Using Approximate_Matching (default)")
endif()

# --------------------------------------------
# Exe
# --------------------------------------------
add_executable(ssm_ged ${SSM_GED_SOURCES})
target_link_libraries(ssm_ged PRIVATE Threads::Threads)
