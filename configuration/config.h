#ifndef _CONFIG_H_
#define _CONFIG_H_

// Shared settings used by both the default AFEC build and paper variants.
#define NDEBUG
#ifndef MATCH_OUTPUT_LIMIT
#define MATCH_OUTPUT_LIMIT 1000000000
#endif

#ifdef SSM_CONFIG_HEADER

#include SSM_CONFIG_HEADER

#else

#define AFEC_ENABLE_SPOKE_FILTERING 1
#define AFEC_ENABLE_BRIDGE_FILTERING 1
#define AFEC_ANCHOR_ORDER_FIXED 0
#define AFEC_ANCHOR_ORDER_DYNAMIC 1
#define AFEC_ANCHOR_ORDER_RANDOM 0
#define AFEC_BW_ALWAYS_BLACK 0
#define AFEC_BW_ALWAYS_WHITE 0
#define AFEC_BW_DYNAMIC 1
#define AFEC_ENABLE_COST_SPLIT 1
#define AFEC_ENABLE_ANCHOR_FRONTIER 1
#define AFEC_USE_FLAT_HASH_MAP 1

#endif // SSM_CONFIG_HEADER

#ifndef AFEC_ENABLE_SPOKE_FILTERING
#define AFEC_ENABLE_SPOKE_FILTERING 1
#endif

#ifndef AFEC_ENABLE_BRIDGE_FILTERING
#define AFEC_ENABLE_BRIDGE_FILTERING 1
#endif

#ifndef AFEC_ANCHOR_ORDER_FIXED
#define AFEC_ANCHOR_ORDER_FIXED 0
#endif

#ifndef AFEC_ANCHOR_ORDER_RANDOM
#define AFEC_ANCHOR_ORDER_RANDOM 0
#endif

#ifndef AFEC_ANCHOR_ORDER_DYNAMIC
#define AFEC_ANCHOR_ORDER_DYNAMIC \
    (!AFEC_ANCHOR_ORDER_FIXED && !AFEC_ANCHOR_ORDER_RANDOM)
#endif

#ifndef AFEC_BW_ALWAYS_BLACK
#define AFEC_BW_ALWAYS_BLACK 0
#endif

#ifndef AFEC_BW_ALWAYS_WHITE
#define AFEC_BW_ALWAYS_WHITE 0
#endif

#ifndef AFEC_BW_DYNAMIC
#define AFEC_BW_DYNAMIC \
    (!AFEC_BW_ALWAYS_BLACK && !AFEC_BW_ALWAYS_WHITE)
#endif

#ifndef AFEC_ENABLE_COST_SPLIT
#define AFEC_ENABLE_COST_SPLIT 1
#endif

#ifndef AFEC_ENABLE_ANCHOR_FRONTIER
#define AFEC_ENABLE_ANCHOR_FRONTIER 1
#endif

#ifndef AFEC_USE_FLAT_HASH_MAP
#define AFEC_USE_FLAT_HASH_MAP 1
#endif

#ifndef MATCH_OUTPUT_LIMIT
#error "MATCH_OUTPUT_LIMIT must be defined by the selected configuration."
#endif

#if AFEC_ENABLE_SPOKE_FILTERING != 0 && AFEC_ENABLE_SPOKE_FILTERING != 1
#error "AFEC_ENABLE_SPOKE_FILTERING must be 0 or 1."
#endif

#if AFEC_ENABLE_BRIDGE_FILTERING != 0 && AFEC_ENABLE_BRIDGE_FILTERING != 1
#error "AFEC_ENABLE_BRIDGE_FILTERING must be 0 or 1."
#endif

#if AFEC_ANCHOR_ORDER_FIXED != 0 && AFEC_ANCHOR_ORDER_FIXED != 1
#error "AFEC_ANCHOR_ORDER_FIXED must be 0 or 1."
#endif

#if AFEC_ANCHOR_ORDER_DYNAMIC != 0 && AFEC_ANCHOR_ORDER_DYNAMIC != 1
#error "AFEC_ANCHOR_ORDER_DYNAMIC must be 0 or 1."
#endif

#if AFEC_ANCHOR_ORDER_RANDOM != 0 && AFEC_ANCHOR_ORDER_RANDOM != 1
#error "AFEC_ANCHOR_ORDER_RANDOM must be 0 or 1."
#endif

#if AFEC_BW_ALWAYS_BLACK != 0 && AFEC_BW_ALWAYS_BLACK != 1
#error "AFEC_BW_ALWAYS_BLACK must be 0 or 1."
#endif

#if AFEC_BW_ALWAYS_WHITE != 0 && AFEC_BW_ALWAYS_WHITE != 1
#error "AFEC_BW_ALWAYS_WHITE must be 0 or 1."
#endif

#if AFEC_BW_DYNAMIC != 0 && AFEC_BW_DYNAMIC != 1
#error "AFEC_BW_DYNAMIC must be 0 or 1."
#endif

#if AFEC_ENABLE_COST_SPLIT != 0 && AFEC_ENABLE_COST_SPLIT != 1
#error "AFEC_ENABLE_COST_SPLIT must be 0 or 1."
#endif

#if AFEC_ENABLE_ANCHOR_FRONTIER != 0 && AFEC_ENABLE_ANCHOR_FRONTIER != 1
#error "AFEC_ENABLE_ANCHOR_FRONTIER must be 0 or 1."
#endif

#if AFEC_USE_FLAT_HASH_MAP != 0 && AFEC_USE_FLAT_HASH_MAP != 1
#error "AFEC_USE_FLAT_HASH_MAP must be 0 or 1."
#endif

#if AFEC_ANCHOR_ORDER_FIXED + AFEC_ANCHOR_ORDER_DYNAMIC + AFEC_ANCHOR_ORDER_RANDOM != 1
#error "Exactly one AFEC ordering strategy must be selected."
#endif

#if AFEC_BW_ALWAYS_BLACK + AFEC_BW_ALWAYS_WHITE + AFEC_BW_DYNAMIC != 1
#error "Exactly one AFEC black--white expansion strategy must be selected."
#endif

#if !AFEC_ENABLE_COST_SPLIT && !AFEC_BW_ALWAYS_BLACK
#error "Disabling Cost-Split requires candidate-wise black expansion (AFE)."
#endif

#if !AFEC_ENABLE_ANCHOR_FRONTIER && AFEC_ENABLE_COST_SPLIT
#error "AFE-NoAF must disable Cost-Split together with anchor-frontier expansion."
#endif

#if MATCH_OUTPUT_LIMIT < 0
#error "MATCH_OUTPUT_LIMIT must be non-negative."
#endif

#endif //_CONFIG_H_
