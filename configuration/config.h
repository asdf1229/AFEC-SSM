#ifndef _CONFIG_H_
#define _CONFIG_H_

// --- cde_edge_ie filtering options ---
#define ENABLE_SPOKE_FILTERING
// #define ENABLE_ONEHOP_FILTERING
#ifndef ONEHOP_FILTER_MISSING_GAP
#define ONEHOP_FILTER_MISSING_GAP 2
#endif
#ifndef ONEHOP_FILTER_MAX_QDEG
#define ONEHOP_FILTER_MAX_QDEG 8
#endif

#if defined(ENABLE_ONEHOP_FILTERING) && !defined(ENABLE_SPOKE_FILTERING)
#error "ENABLE_ONEHOP_FILTERING requires ENABLE_SPOKE_FILTERING because OneHop consumes spoke-stage records."
#endif

// --- cde_edge_ie lower-bound options ---
// Disabled by default: current experiments show the search-time lower bound
// costs more than it prunes.
// #define CDE_LB_LIGHTWEIGHT_SPOKE

// --- other options ---
#define NDEBUG

// --- output options ---
// 0 means unlimited. A positive value stops recursive enumeration after that
// many results have been found.
#ifndef MATCH_OUTPUT_LIMIT
#define MATCH_OUTPUT_LIMIT 0
#endif
#if MATCH_OUTPUT_LIMIT < 0
#error "MATCH_OUTPUT_LIMIT must be non-negative."
#endif

// --- cde_edge_ie anchor-support options ---
// Cached anchor-support is the default. It reuses support across child states
// and only recomputes entries marked dirty by search-state changes.
#define CDE_EDGE_IE_CACHE_ANCHOR_SUPPORT
// Recompute anchor-support from the current search state when it is scored.
// This is the exact baseline for comparing cache behavior.
// #define CDE_EDGE_IE_RECOMPUTE_ANCHOR_SUPPORT
// #define ENABLE_EXCLUDED_EDGE_SUPPORT

#endif //_CONFIG_H_
