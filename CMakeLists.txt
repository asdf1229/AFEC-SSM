cmake_minimum_required(VERSION 3.10)

project(SSM_GED LANGUAGES CXX)

add_compile_options(-O3 -g -Wall)

include_directories(
    ${PROJECT_SOURCE_DIR}
)

find_package(Threads REQUIRED)

set(SSM_GED_COMMON_SOURCES
    ${PROJECT_SOURCE_DIR}/graph/graph.cpp
    ${PROJECT_SOURCE_DIR}/matching/main.cpp
    ${PROJECT_SOURCE_DIR}/matching/run_matching.cpp
)

file(GLOB SSM_GED_ALGORITHM_SOURCES
    ${PROJECT_SOURCE_DIR}/matching/algorithms/*.cpp
)
list(SORT SSM_GED_ALGORITHM_SOURCES)

if (NOT SSM_GED_ALGORITHM_SOURCES)
    message(FATAL_ERROR "No algorithm entry files found under matching/algorithms.")
endif()

set(SSM_GED_EXECUTABLES "")

foreach(ALGO_SOURCE ${SSM_GED_ALGORITHM_SOURCES})
    get_filename_component(ALGO_NAME ${ALGO_SOURCE} NAME_WE)
    set(TARGET_NAME "ssm_ged_${ALGO_NAME}")

    add_executable(${TARGET_NAME}
        ${SSM_GED_COMMON_SOURCES}
        ${ALGO_SOURCE}
    )
    target_link_libraries(${TARGET_NAME} PRIVATE Threads::Threads)

    list(APPEND SSM_GED_EXECUTABLES ${TARGET_NAME})
endforeach()

message(STATUS "Configured SSM-GED executables: ${SSM_GED_EXECUTABLES}")
