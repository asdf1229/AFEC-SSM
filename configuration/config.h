#ifndef _CONFIG_H_
#define _CONFIG_H_

// Shared settings used by both default and ablation configurations.
#define NDEBUG
#ifndef MATCH_OUTPUT_LIMIT
#define MATCH_OUTPUT_LIMIT 1000000000
#endif

#ifdef SSM_CONFIG_HEADER

#include SSM_CONFIG_HEADER

#else

#define AFEE_ENABLE_SPOKE_FILTERING 1
#define AFEE_ENABLE_BRIDGE_FILTERING 1
#define AFEE_FIXED_ORDER 0
#define AFEE_DYNAMIC_ORDER 1
#define AFEE_RANDOM_ORDER 0
#define AFEE_COLOR_ALL_BLACK 0
#define AFEE_COLOR_ALL_WHITE 0
#define AFEE_DYNAMIC_COLOR 1
#define AFEE_ENABLE_SPLIT 1
#define AFEE_USE_FLAT_HASH_MAP 1

#endif // SSM_CONFIG_HEADER

#ifndef AFEE_ENABLE_SPOKE_FILTERING
#define AFEE_ENABLE_SPOKE_FILTERING 1
#endif

#ifndef AFEE_ENABLE_BRIDGE_FILTERING
#define AFEE_ENABLE_BRIDGE_FILTERING 1
#endif

#ifndef AFEE_FIXED_ORDER
#define AFEE_FIXED_ORDER 0
#endif

#ifndef AFEE_RANDOM_ORDER
#define AFEE_RANDOM_ORDER 0
#endif

#ifndef AFEE_DYNAMIC_ORDER
#define AFEE_DYNAMIC_ORDER \
    (!AFEE_FIXED_ORDER && !AFEE_RANDOM_ORDER)
#endif

#ifndef AFEE_COLOR_ALL_BLACK
#define AFEE_COLOR_ALL_BLACK 0
#endif

#ifndef AFEE_COLOR_ALL_WHITE
#define AFEE_COLOR_ALL_WHITE 0
#endif

#ifndef AFEE_DYNAMIC_COLOR
#define AFEE_DYNAMIC_COLOR \
    (!AFEE_COLOR_ALL_BLACK && !AFEE_COLOR_ALL_WHITE)
#endif

#ifndef AFEE_ENABLE_SPLIT
#define AFEE_ENABLE_SPLIT 1
#endif

#ifndef AFEE_USE_FLAT_HASH_MAP
#define AFEE_USE_FLAT_HASH_MAP 1
#endif

#ifndef MATCH_OUTPUT_LIMIT
#error "MATCH_OUTPUT_LIMIT must be defined by the selected configuration."
#endif

#if AFEE_ENABLE_SPOKE_FILTERING != 0 && AFEE_ENABLE_SPOKE_FILTERING != 1
#error "AFEE_ENABLE_SPOKE_FILTERING must be 0 or 1."
#endif

#if AFEE_ENABLE_BRIDGE_FILTERING != 0 && AFEE_ENABLE_BRIDGE_FILTERING != 1
#error "AFEE_ENABLE_BRIDGE_FILTERING must be 0 or 1."
#endif

#if AFEE_FIXED_ORDER != 0 && AFEE_FIXED_ORDER != 1
#error "AFEE_FIXED_ORDER must be 0 or 1."
#endif

#if AFEE_DYNAMIC_ORDER != 0 && AFEE_DYNAMIC_ORDER != 1
#error "AFEE_DYNAMIC_ORDER must be 0 or 1."
#endif

#if AFEE_RANDOM_ORDER != 0 && AFEE_RANDOM_ORDER != 1
#error "AFEE_RANDOM_ORDER must be 0 or 1."
#endif

#if AFEE_COLOR_ALL_BLACK != 0 && AFEE_COLOR_ALL_BLACK != 1
#error "AFEE_COLOR_ALL_BLACK must be 0 or 1."
#endif

#if AFEE_COLOR_ALL_WHITE != 0 && AFEE_COLOR_ALL_WHITE != 1
#error "AFEE_COLOR_ALL_WHITE must be 0 or 1."
#endif

#if AFEE_DYNAMIC_COLOR != 0 && AFEE_DYNAMIC_COLOR != 1
#error "AFEE_DYNAMIC_COLOR must be 0 or 1."
#endif

#if AFEE_ENABLE_SPLIT != 0 && AFEE_ENABLE_SPLIT != 1
#error "AFEE_ENABLE_SPLIT must be 0 or 1."
#endif

#if AFEE_USE_FLAT_HASH_MAP != 0 && AFEE_USE_FLAT_HASH_MAP != 1
#error "AFEE_USE_FLAT_HASH_MAP must be 0 or 1."
#endif

#if AFEE_FIXED_ORDER + AFEE_DYNAMIC_ORDER + AFEE_RANDOM_ORDER != 1
#error "Exactly one AFEE ordering strategy must be selected."
#endif

#if AFEE_COLOR_ALL_BLACK + AFEE_COLOR_ALL_WHITE + AFEE_DYNAMIC_COLOR != 1
#error "Exactly one AFEE coloring strategy must be selected."
#endif

#if MATCH_OUTPUT_LIMIT < 0
#error "MATCH_OUTPUT_LIMIT must be non-negative."
#endif

#endif //_CONFIG_H_
