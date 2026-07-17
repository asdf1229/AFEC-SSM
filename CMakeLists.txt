cmake_minimum_required(VERSION 3.16)

project(SSM LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# CMake adds -DNDEBUG to its non-Debug build-type flags by default.  NDEBUG is
# intentionally owned by configuration/config.h in this project, so remove it
# from every built-in optimization configuration.
foreach(SSM_BUILD_TYPE RELEASE RELWITHDEBINFO MINSIZEREL)
    string(REPLACE "-DNDEBUG" "" SSM_CXX_FLAGS
        "${CMAKE_CXX_FLAGS_${SSM_BUILD_TYPE}}")
    set(CMAKE_CXX_FLAGS_${SSM_BUILD_TYPE} "${SSM_CXX_FLAGS}" CACHE STRING
        "Flags used by the C++ compiler during ${SSM_BUILD_TYPE} builds." FORCE)
endforeach()

add_compile_options(-O3 -g -Wall)

set(SSM_COMPILE_DEFINITIONS "" CACHE STRING
    "Extra semicolon-separated compile definitions for SSM-GED targets.")
set(SSM_CONFIG_HEADER "" CACHE STRING
    "Optional configuration header selected instead of configuration/config.h defaults.")

include_directories(
    ${PROJECT_SOURCE_DIR}
)

find_package(Threads REQUIRED)
find_package(absl REQUIRED CONFIG)

message(STATUS "C++ compiler: ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Abseil package: ${absl_DIR}")

set(SSM_COMMON_SOURCES
    ${PROJECT_SOURCE_DIR}/graph/graph.cpp
    ${PROJECT_SOURCE_DIR}/matching/main.cpp
    ${PROJECT_SOURCE_DIR}/matching/run_matching.cpp
)

file(GLOB SSM_ALGORITHM_SOURCES CONFIGURE_DEPENDS
    ${PROJECT_SOURCE_DIR}/matching/algorithms/*.cpp
)
list(APPEND SSM_ALGORITHM_SOURCES
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afee/afee.cpp
)
list(SORT SSM_ALGORITHM_SOURCES)

set(SSM_AFEE_SOURCES
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afee/branching.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afee/buffers.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afee/candidate_index.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afee/filtering.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afee/ordering.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afee/state_undo.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afee/terminal_tail.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afee/white_candidates.cpp
)

if (NOT SSM_ALGORITHM_SOURCES)
    message(FATAL_ERROR "No algorithm entry files found under matching/algorithms.")
endif()

set(SSM_EXECUTABLES "")

macro(add_ssm_executable TARGET_NAME ALGO_SOURCE CONFIG_HEADER)
    add_executable(${TARGET_NAME}
        ${SSM_COMMON_SOURCES}
        ${ALGO_SOURCE}
        ${ARGN}
    )
    if (SSM_COMPILE_DEFINITIONS)
        target_compile_definitions(${TARGET_NAME} PRIVATE ${SSM_COMPILE_DEFINITIONS})
    endif()
    if (NOT "${CONFIG_HEADER}" STREQUAL "")
        target_compile_definitions(${TARGET_NAME} PRIVATE
            "SSM_CONFIG_HEADER=\"${CONFIG_HEADER}\"")
    endif()
    target_link_libraries(${TARGET_NAME} PRIVATE Threads::Threads absl::hash absl::raw_hash_set absl::base)
    list(APPEND SSM_EXECUTABLES ${TARGET_NAME})
endmacro()

macro(set_afee_identity TARGET_NAME CONFIG_NAME)
    if ("${CONFIG_NAME}" STREQUAL "afee")
        set(AFEE_ALGORITHM_KEY "afee")
        set(AFEE_DISPLAY_NAME "AFEE")
    else()
        set(AFEE_ALGORITHM_KEY "afee_${CONFIG_NAME}")
        string(TOUPPER "${CONFIG_NAME}" AFEE_CONFIG_DISPLAY_NAME)
        set(AFEE_DISPLAY_NAME "AFEE_${AFEE_CONFIG_DISPLAY_NAME}")
    endif()
    target_compile_definitions(${TARGET_NAME} PRIVATE
        "AFEE_ALGORITHM_KEY=\"${AFEE_ALGORITHM_KEY}\""
        "AFEE_DISPLAY_NAME=\"${AFEE_DISPLAY_NAME}\""
    )
endmacro()

macro(add_afee_executable ALGO_SOURCE CONFIG_HEADER)
    get_filename_component(CONFIG_NAME "${CONFIG_HEADER}" NAME_WE)
    if (CONFIG_NAME STREQUAL "afee")
        set(AFEE_TARGET_NAME "ssm_afee")
    else()
        set(AFEE_TARGET_NAME "ssm_afee_${CONFIG_NAME}")
    endif()
    add_ssm_executable(
        "${AFEE_TARGET_NAME}"
        "${ALGO_SOURCE}"
        "${CONFIG_HEADER}"
        ${SSM_AFEE_SOURCES}
    )
    set_afee_identity("${AFEE_TARGET_NAME}" "${CONFIG_NAME}")
endmacro()

foreach(ALGO_SOURCE ${SSM_ALGORITHM_SOURCES})
    get_filename_component(ALGO_FILE_NAME ${ALGO_SOURCE} NAME_WE)
    get_filename_component(ALGO_DIR ${ALGO_SOURCE} DIRECTORY)
    get_filename_component(ALGO_DIR_NAME ${ALGO_DIR} NAME)
    if (ALGO_FILE_NAME STREQUAL ALGO_DIR_NAME)
        set(ALGO_NAME ${ALGO_DIR_NAME})
    else()
        set(ALGO_NAME ${ALGO_FILE_NAME})
    endif()

    if (ALGO_NAME STREQUAL "afee")
        if (SSM_CONFIG_HEADER)
            add_afee_executable("${ALGO_SOURCE}" "${SSM_CONFIG_HEADER}")
        else()
            file(GLOB AFEE_ABLATION_CONFIGS CONFIGURE_DEPENDS
                RELATIVE "${PROJECT_SOURCE_DIR}"
                "${PROJECT_SOURCE_DIR}/configuration/afee_ablation_configs/*.h"
            )
            list(SORT AFEE_ABLATION_CONFIGS)
            if (NOT AFEE_ABLATION_CONFIGS)
                message(FATAL_ERROR
                    "No AFEE ablation configs found under configuration/afee_ablation_configs.")
            endif()

            foreach(CONFIG_HEADER ${AFEE_ABLATION_CONFIGS})
                add_afee_executable("${ALGO_SOURCE}" "${CONFIG_HEADER}")
            endforeach()
        endif()
    elseif(ALGO_NAME STREQUAL "s3and")
        add_ssm_executable(
            "ssm_s3and_ssm"
            "${ALGO_SOURCE}"
            ""
        )
    else()
        add_ssm_executable(
            "ssm_${ALGO_NAME}"
            "${ALGO_SOURCE}"
            ""
        )
    endif()
endforeach()

message(STATUS "Configured SSM-GED executables: ${SSM_EXECUTABLES}")
