#ifndef _CONFIG_H_
#define _CONFIG_H_

// --- cde_edge_ie filtering options ---
#ifndef DISABLE_SPOKE_FILTERING
#define ENABLE_SPOKE_FILTERING
#endif
// #define ENABLE_ONEHOP_FILTERING
#ifdef DISABLE_ONEHOP_FILTERING
#undef ENABLE_ONEHOP_FILTERING
#endif
#ifndef ONEHOP_FILTER_MISSING_GAP
#define ONEHOP_FILTER_MISSING_GAP 2
#endif
#ifndef ONEHOP_FILTER_MAX_QDEG
#define ONEHOP_FILTER_MAX_QDEG 8
#endif

#if defined(ENABLE_ONEHOP_FILTERING) && !defined(ENABLE_SPOKE_FILTERING)
#error "ENABLE_ONEHOP_FILTERING requires ENABLE_SPOKE_FILTERING because OneHop consumes spoke-stage records."
#endif

// --- other options ---
#define NDEBUG

// --- output options ---
// 0 means unlimited. A positive value stops recursive enumeration after that
// many results have been found.
#ifndef MATCH_OUTPUT_LIMIT
#define MATCH_OUTPUT_LIMIT 1000
#endif
#if MATCH_OUTPUT_LIMIT < 0
#error "MATCH_OUTPUT_LIMIT must be non-negative."
#endif

// #define ENABLE_EXCLUDED_EDGE_SUPPORT
// #define ENABLE_CAND_STATS

#endif //_CONFIG_H_
