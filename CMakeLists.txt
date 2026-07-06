cmake_minimum_required(VERSION 3.10)

project(SSM_GED LANGUAGES CXX)

add_compile_options(-O3 -g -Wall)

set(SSM_GED_COMPILE_DEFINITIONS "" CACHE STRING
    "Extra semicolon-separated compile definitions for SSM-GED targets.")
set(SSM_GED_CONFIG_HEADER "" CACHE STRING
    "Optional configuration header selected instead of configuration/config.h defaults.")

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
list(APPEND SSM_GED_ALGORITHM_SOURCES
    ${PROJECT_SOURCE_DIR}/matching/algorithms/cde_black_white/cde_black_white.cpp
)
list(SORT SSM_GED_ALGORITHM_SOURCES)

set(SSM_GED_CDE_BLACK_WHITE_SOURCES
    ${PROJECT_SOURCE_DIR}/matching/algorithms/cde_black_white/branching.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/cde_black_white/buffers.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/cde_black_white/candidate_index.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/cde_black_white/filtering.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/cde_black_white/ordering.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/cde_black_white/state_undo.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/cde_black_white/terminal_tail.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/cde_black_white/white_candidates.cpp
)

if (NOT SSM_GED_ALGORITHM_SOURCES)
    message(FATAL_ERROR "No algorithm entry files found under matching/algorithms.")
endif()

set(SSM_GED_EXECUTABLES "")

macro(add_ssm_ged_executable TARGET_NAME ALGO_SOURCE CONFIG_HEADER)
    add_executable(${TARGET_NAME}
        ${SSM_GED_COMMON_SOURCES}
        ${ALGO_SOURCE}
        ${ARGN}
    )
    if (SSM_GED_COMPILE_DEFINITIONS)
        target_compile_definitions(${TARGET_NAME} PRIVATE ${SSM_GED_COMPILE_DEFINITIONS})
    endif()
    if (NOT "${CONFIG_HEADER}" STREQUAL "")
        target_compile_definitions(${TARGET_NAME} PRIVATE
            "SSM_GED_CONFIG_HEADER=\"${CONFIG_HEADER}\"")
    endif()
    target_link_libraries(${TARGET_NAME} PRIVATE Threads::Threads)
    list(APPEND SSM_GED_EXECUTABLES ${TARGET_NAME})
endmacro()

macro(set_cde_edge_ie_identity TARGET_NAME CONFIG_NAME)
    if ("${CONFIG_NAME}" STREQUAL "base")
        set(CDE_ALGORITHM_KEY "cde_edge_ie")
        set(CDE_DISPLAY_NAME "CDE-Edge-IE")
    else()
        set(CDE_ALGORITHM_KEY "cde_edge_ie_${CONFIG_NAME}")
        set(CDE_DISPLAY_NAME "CDE-Edge-IE (${CONFIG_NAME})")
    endif()
    target_compile_definitions(${TARGET_NAME} PRIVATE
        "CDE_EDGE_IE_ALGORITHM_KEY=\"${CDE_ALGORITHM_KEY}\""
        "CDE_EDGE_IE_DISPLAY_NAME=\"${CDE_DISPLAY_NAME}\""
    )
endmacro()

macro(set_cde_black_white_identity TARGET_NAME CONFIG_NAME)
    set(CDE_BLACK_WHITE_ALGORITHM_KEY "cde_black_white_${CONFIG_NAME}")
    set(CDE_BLACK_WHITE_DISPLAY_NAME "CDE-Black-White (${CONFIG_NAME})")
    target_compile_definitions(${TARGET_NAME} PRIVATE
        "CDE_BLACK_WHITE_ALGORITHM_KEY=\"${CDE_BLACK_WHITE_ALGORITHM_KEY}\""
        "CDE_BLACK_WHITE_DISPLAY_NAME=\"${CDE_BLACK_WHITE_DISPLAY_NAME}\""
    )
endmacro()

foreach(ALGO_SOURCE ${SSM_GED_ALGORITHM_SOURCES})
    get_filename_component(ALGO_FILE_NAME ${ALGO_SOURCE} NAME_WE)
    get_filename_component(ALGO_DIR ${ALGO_SOURCE} DIRECTORY)
    get_filename_component(ALGO_DIR_NAME ${ALGO_DIR} NAME)
    if (ALGO_FILE_NAME STREQUAL ALGO_DIR_NAME)
        set(ALGO_NAME ${ALGO_DIR_NAME})
    else()
        set(ALGO_NAME ${ALGO_FILE_NAME})
    endif()

    if (ALGO_NAME STREQUAL "cde_edge_ie")
        if (SSM_GED_CONFIG_HEADER)
            # Explicit single-config mode is used by scripts/run_cde_ablations.py.
            get_filename_component(CONFIG_NAME "${SSM_GED_CONFIG_HEADER}" NAME_WE)
            add_ssm_ged_executable(
                "ssm_ged_cde_edge_ie"
                "${ALGO_SOURCE}"
                "${SSM_GED_CONFIG_HEADER}"
            )
            set_cde_edge_ie_identity("ssm_ged_cde_edge_ie" "${CONFIG_NAME}")
        else()
            # Normal build.sh mode: every config header becomes one executable.
            file(GLOB CDE_ABLATION_CONFIGS CONFIGURE_DEPENDS
                RELATIVE "${PROJECT_SOURCE_DIR}"
                "${PROJECT_SOURCE_DIR}/configuration/cde_ablation_configs/*.h"
            )
            list(SORT CDE_ABLATION_CONFIGS)
            if (NOT CDE_ABLATION_CONFIGS)
                message(FATAL_ERROR
                    "No CDE ablation configs found under configuration/cde_ablation_configs.")
            endif()

            set(CDE_BASE_CONFIG_FOUND FALSE)
            foreach(CONFIG_HEADER ${CDE_ABLATION_CONFIGS})
                get_filename_component(CONFIG_NAME "${CONFIG_HEADER}" NAME_WE)
                if (CONFIG_NAME STREQUAL "base")
                    set(TARGET_NAME "ssm_ged_cde_edge_ie")
                    set(CDE_BASE_CONFIG_FOUND TRUE)
                else()
                    set(TARGET_NAME "ssm_ged_cde_edge_ie_${CONFIG_NAME}")
                endif()
                add_ssm_ged_executable(
                    "${TARGET_NAME}"
                    "${ALGO_SOURCE}"
                    "${CONFIG_HEADER}"
                )
                set_cde_edge_ie_identity("${TARGET_NAME}" "${CONFIG_NAME}")
            endforeach()

            if (NOT CDE_BASE_CONFIG_FOUND)
                message(FATAL_ERROR
                    "CDE ablation configs must include configuration/cde_ablation_configs/base.h.")
            endif()
        endif()
    elseif (ALGO_NAME STREQUAL "cde_black_white")
        if (SSM_GED_CONFIG_HEADER)
            get_filename_component(CONFIG_NAME "${SSM_GED_CONFIG_HEADER}" NAME_WE)
            add_ssm_ged_executable(
                "ssm_ged_cde_black_white_${CONFIG_NAME}"
                "${ALGO_SOURCE}"
                "${SSM_GED_CONFIG_HEADER}"
                ${SSM_GED_CDE_BLACK_WHITE_SOURCES}
            )
            set_cde_black_white_identity("ssm_ged_cde_black_white_${CONFIG_NAME}" "${CONFIG_NAME}")
        else()
            file(GLOB CDE_BLACK_WHITE_ABLATION_CONFIGS CONFIGURE_DEPENDS
                RELATIVE "${PROJECT_SOURCE_DIR}"
                "${PROJECT_SOURCE_DIR}/configuration/cde_black_white_ablation_configs/*.h"
            )
            list(SORT CDE_BLACK_WHITE_ABLATION_CONFIGS)
            if (NOT CDE_BLACK_WHITE_ABLATION_CONFIGS)
                message(FATAL_ERROR
                    "No CDE black-white ablation configs found under configuration/cde_black_white_ablation_configs.")
            endif()

            foreach(CONFIG_HEADER ${CDE_BLACK_WHITE_ABLATION_CONFIGS})
                get_filename_component(CONFIG_NAME "${CONFIG_HEADER}" NAME_WE)
                set(TARGET_NAME "ssm_ged_cde_black_white_${CONFIG_NAME}")
                add_ssm_ged_executable(
                    "${TARGET_NAME}"
                    "${ALGO_SOURCE}"
                    "${CONFIG_HEADER}"
                    ${SSM_GED_CDE_BLACK_WHITE_SOURCES}
                )
                set_cde_black_white_identity("${TARGET_NAME}" "${CONFIG_NAME}")
            endforeach()
        endif()
    else()
        add_ssm_ged_executable(
            "ssm_ged_${ALGO_NAME}"
            "${ALGO_SOURCE}"
            ""
        )
    endif()
endforeach()

message(STATUS "Configured SSM-GED executables: ${SSM_GED_EXECUTABLES}")
