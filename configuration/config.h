#ifndef _CONFIG_H_
#define _CONFIG_H_

// Shared settings used by both default and ablation configurations.
#define NDEBUG
#define MATCH_OUTPUT_LIMIT 1000000
#define CDE_EDGE_IE_TOPK_SUPPORT_DECAY_GAMMA 0.9

#ifdef SSM_GED_CONFIG_HEADER

#include SSM_GED_CONFIG_HEADER

#else

// #define ENABLE_EXCLUDED_EDGE_SUPPORT
// #define ENABLE_CAND_STATS

#define CDE_EDGE_IE_TERMINAL_BUCKETS_DEFAULT 0
#define CDE_EDGE_IE_ENABLE_SPOKE_FILTERING 1
#define CDE_EDGE_IE_ENABLE_BRIDGE_FILTERING 1
#define CDE_EDGE_IE_FIXED_ORDER 0
#define CDE_EDGE_IE_TOPK_SUPPORT_DECAY 0

#endif // SSM_GED_CONFIG_HEADER

#ifndef MATCH_OUTPUT_LIMIT
#error "MATCH_OUTPUT_LIMIT must be defined by the selected configuration."
#endif

#ifndef CDE_EDGE_IE_TERMINAL_BUCKETS_DEFAULT
#error "CDE_EDGE_IE_TERMINAL_BUCKETS_DEFAULT must be defined by the selected configuration."
#endif

#ifndef CDE_EDGE_IE_ENABLE_SPOKE_FILTERING
#error "CDE_EDGE_IE_ENABLE_SPOKE_FILTERING must be defined by the selected configuration."
#endif

#ifndef CDE_EDGE_IE_ENABLE_BRIDGE_FILTERING
#error "CDE_EDGE_IE_ENABLE_BRIDGE_FILTERING must be defined by the selected configuration."
#endif

#ifndef CDE_EDGE_IE_FIXED_ORDER
#error "CDE_EDGE_IE_FIXED_ORDER must be defined by the selected configuration."
#endif

#ifndef CDE_EDGE_IE_TOPK_SUPPORT_DECAY
#error "CDE_EDGE_IE_TOPK_SUPPORT_DECAY must be defined by the selected configuration."
#endif

#ifndef CDE_EDGE_IE_TOPK_SUPPORT_DECAY_GAMMA
#error "CDE_EDGE_IE_TOPK_SUPPORT_DECAY_GAMMA must be defined by the selected configuration."
#endif

#if CDE_EDGE_IE_TOPK_SUPPORT_DECAY && CDE_EDGE_IE_FIXED_ORDER
#error "CDE_EDGE_IE_TOPK_SUPPORT_DECAY and CDE_EDGE_IE_FIXED_ORDER are mutually exclusive."
#endif

#if CDE_EDGE_IE_TERMINAL_BUCKETS_DEFAULT != 0 && CDE_EDGE_IE_TERMINAL_BUCKETS_DEFAULT != 1
#error "CDE_EDGE_IE_TERMINAL_BUCKETS_DEFAULT must be 0 or 1."
#endif

#if CDE_EDGE_IE_ENABLE_SPOKE_FILTERING != 0 && CDE_EDGE_IE_ENABLE_SPOKE_FILTERING != 1
#error "CDE_EDGE_IE_ENABLE_SPOKE_FILTERING must be 0 or 1."
#endif

#if CDE_EDGE_IE_ENABLE_BRIDGE_FILTERING != 0 && CDE_EDGE_IE_ENABLE_BRIDGE_FILTERING != 1
#error "CDE_EDGE_IE_ENABLE_BRIDGE_FILTERING must be 0 or 1."
#endif

#if CDE_EDGE_IE_FIXED_ORDER != 0 && CDE_EDGE_IE_FIXED_ORDER != 1
#error "CDE_EDGE_IE_FIXED_ORDER must be 0 or 1."
#endif

#if CDE_EDGE_IE_TOPK_SUPPORT_DECAY != 0 && CDE_EDGE_IE_TOPK_SUPPORT_DECAY != 1
#error "CDE_EDGE_IE_TOPK_SUPPORT_DECAY must be 0 or 1."
#endif

#if MATCH_OUTPUT_LIMIT < 0
#error "MATCH_OUTPUT_LIMIT must be non-negative."
#endif

#endif //_CONFIG_H_
