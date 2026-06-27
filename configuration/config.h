#ifndef _CONFIG_H_
#define _CONFIG_H_

#ifdef SSM_GED_CONFIG_HEADER

#include SSM_GED_CONFIG_HEADER

#else

// --- cde_edge_ie filtering options ---
#ifndef DISABLE_SPOKE_FILTERING
#define ENABLE_SPOKE_FILTERING
#endif
// #define ENABLE_ONEHOP_FILTERING
#ifdef DISABLE_ONEHOP_FILTERING
#undef ENABLE_ONEHOP_FILTERING
#endif

// --- other options ---
#define NDEBUG

// --- output options ---
// 0 means unlimited. A positive value stops recursive enumeration after that
// many results have been found.

// #define ENABLE_EXCLUDED_EDGE_SUPPORT
// #define ENABLE_CAND_STATS

#define CDE_EDGE_IE_TERMINAL_BUCKETS_DEFAULT 1
#define CDE_EDGE_IE_ENABLE_BRIDGE_FILTERING 1
#define CDE_EDGE_IE_FIXED_ORDER 0
#define CDE_EDGE_IE_TOPK_SUPPORT_DECAY 0
#define CDE_EDGE_IE_TOPK_SUPPORT_DECAY_GAMMA 0.9

#endif // SSM_GED_CONFIG_HEADER

#ifndef ONEHOP_FILTER_MISSING_GAP
#define ONEHOP_FILTER_MISSING_GAP 2
#endif

#ifndef ONEHOP_FILTER_MAX_QDEG
#define ONEHOP_FILTER_MAX_QDEG 8
#endif

#ifndef MATCH_OUTPUT_LIMIT
#define MATCH_OUTPUT_LIMIT 1000
#endif

#ifndef CDE_EDGE_IE_TERMINAL_BUCKETS_DEFAULT
#define CDE_EDGE_IE_TERMINAL_BUCKETS_DEFAULT 1
#endif

#ifndef CDE_EDGE_IE_ENABLE_BRIDGE_FILTERING
#define CDE_EDGE_IE_ENABLE_BRIDGE_FILTERING 1
#endif

#ifndef CDE_EDGE_IE_FIXED_ORDER
#define CDE_EDGE_IE_FIXED_ORDER 0
#endif

#ifndef CDE_EDGE_IE_TOPK_SUPPORT_DECAY
#define CDE_EDGE_IE_TOPK_SUPPORT_DECAY 0
#endif

#ifndef CDE_EDGE_IE_TOPK_SUPPORT_DECAY_GAMMA
#define CDE_EDGE_IE_TOPK_SUPPORT_DECAY_GAMMA 0.9
#endif

#if CDE_EDGE_IE_TOPK_SUPPORT_DECAY && CDE_EDGE_IE_FIXED_ORDER
#error "CDE_EDGE_IE_TOPK_SUPPORT_DECAY and CDE_EDGE_IE_FIXED_ORDER are mutually exclusive."
#endif

#if defined(ENABLE_ONEHOP_FILTERING) && !defined(ENABLE_SPOKE_FILTERING)
#error "ENABLE_ONEHOP_FILTERING requires ENABLE_SPOKE_FILTERING because OneHop consumes spoke-stage records."
#endif

#if MATCH_OUTPUT_LIMIT < 0
#error "MATCH_OUTPUT_LIMIT must be non-negative."
#endif

#endif //_CONFIG_H_
