cmake_minimum_required(VERSION 3.16)

project(AFEC LANGUAGES CXX)

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
    "Extra semicolon-separated compile definitions for AFEC suite targets.")
set(SSM_CONFIG_HEADER "" CACHE STRING
    "Optional AFEC variant header; when set, only that AFEC variant is built.")

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

file(GLOB SSM_BASELINE_SOURCES CONFIGURE_DEPENDS
    ${PROJECT_SOURCE_DIR}/matching/algorithms/*.cpp
)
list(SORT SSM_BASELINE_SOURCES)

set(AFEC_ENTRY_SOURCE
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afec/afec.cpp
)

set(AFEC_SOURCES
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afec/branching.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afec/buffers.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afec/candidate_index.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afec/filtering.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afec/no_anchor_frontier.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afec/ordering.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afec/state_undo.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afec/terminal_tail.cpp
    ${PROJECT_SOURCE_DIR}/matching/algorithms/afec/white_candidates.cpp
)

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
    target_link_libraries(${TARGET_NAME} PRIVATE
        Threads::Threads absl::hash absl::raw_hash_set absl::base)
    list(APPEND SSM_EXECUTABLES ${TARGET_NAME})
endmacro()

macro(resolve_afec_identity CONFIG_HEADER)
    get_filename_component(AFEC_CONFIG_NAME "${CONFIG_HEADER}" NAME_WE)
    if (AFEC_CONFIG_NAME STREQUAL "afec")
        set(AFEC_ALGORITHM_KEY "afec")
        set(AFEC_DISPLAY_NAME "AFEC")
    elseif(AFEC_CONFIG_NAME STREQUAL "afe")
        set(AFEC_ALGORITHM_KEY "afe")
        set(AFEC_DISPLAY_NAME "AFE")
    elseif(AFEC_CONFIG_NAME STREQUAL "afe_no_af")
        set(AFEC_ALGORITHM_KEY "afe_no_af")
        set(AFEC_DISPLAY_NAME "AFE-NoAF")
    elseif(AFEC_CONFIG_NAME STREQUAL "afec_white")
        set(AFEC_ALGORITHM_KEY "afec_white")
        set(AFEC_DISPLAY_NAME "AFEC-White")
    elseif(AFEC_CONFIG_NAME STREQUAL "afec_black")
        set(AFEC_ALGORITHM_KEY "afec_black")
        set(AFEC_DISPLAY_NAME "AFEC-Black")
    elseif(AFEC_CONFIG_NAME STREQUAL "afec_fixed")
        set(AFEC_ALGORITHM_KEY "afec_fixed")
        set(AFEC_DISPLAY_NAME "AFEC-Fixed")
    elseif(AFEC_CONFIG_NAME STREQUAL "afec_random")
        set(AFEC_ALGORITHM_KEY "afec_random")
        set(AFEC_DISPLAY_NAME "AFEC-Random")
    else()
        set(AFEC_ALGORITHM_KEY "afec_${AFEC_CONFIG_NAME}")
        string(REPLACE "_" "-" AFEC_DISPLAY_SUFFIX "${AFEC_CONFIG_NAME}")
        set(AFEC_DISPLAY_NAME "AFEC-${AFEC_DISPLAY_SUFFIX}")
    endif()
endmacro()

macro(add_afec_variant CONFIG_HEADER ALGORITHM_KEY DISPLAY_NAME)
    set(AFEC_TARGET_NAME "ssm_${ALGORITHM_KEY}")
    add_ssm_executable(
        "${AFEC_TARGET_NAME}"
        "${AFEC_ENTRY_SOURCE}"
        "${CONFIG_HEADER}"
        ${AFEC_SOURCES}
    )
    target_compile_definitions(${AFEC_TARGET_NAME} PRIVATE
        "AFEC_ALGORITHM_KEY=\"${ALGORITHM_KEY}\""
        "AFEC_DISPLAY_NAME=\"${DISPLAY_NAME}\""
    )
endmacro()

if (SSM_CONFIG_HEADER)
    resolve_afec_identity("${SSM_CONFIG_HEADER}")
    add_afec_variant(
        "${SSM_CONFIG_HEADER}" "${AFEC_ALGORITHM_KEY}" "${AFEC_DISPLAY_NAME}")
else()
    # Keep the paper's seven public configurations explicit. This prevents
    # unrelated diagnostic headers from silently becoming experiment variants.
    add_afec_variant("configuration/afec_variants/afec.h" "afec" "AFEC")
    add_afec_variant("configuration/afec_variants/afe.h" "afe" "AFE")
    add_afec_variant(
        "configuration/afec_variants/afe_no_af.h" "afe_no_af" "AFE-NoAF")
    add_afec_variant(
        "configuration/afec_variants/afec_white.h" "afec_white" "AFEC-White")
    add_afec_variant(
        "configuration/afec_variants/afec_black.h" "afec_black" "AFEC-Black")
    add_afec_variant(
        "configuration/afec_variants/afec_fixed.h" "afec_fixed" "AFEC-Fixed")
    add_afec_variant(
        "configuration/afec_variants/afec_random.h" "afec_random" "AFEC-Random")
endif()

foreach(ALGO_SOURCE ${SSM_BASELINE_SOURCES})
    get_filename_component(ALGO_FILE_NAME ${ALGO_SOURCE} NAME_WE)
    set(ALGO_NAME ${ALGO_FILE_NAME})

    if(ALGO_NAME STREQUAL "s3and")
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

message(STATUS "Configured AFEC suite executables: ${SSM_EXECUTABLES}")
