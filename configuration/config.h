#ifndef _CONFIG_H_
#define _CONFIG_H_

// --- cde_edge_ie filtering options ---
#define ENABLE_SPOKE_FILTERING
#define ENABLE_ONEHOP_FILTERING
#ifndef ONEHOP_FILTER_MISSING_GAP
#define ONEHOP_FILTER_MISSING_GAP 2
#endif

#if defined(ENABLE_ONEHOP_FILTERING) && !defined(ENABLE_SPOKE_FILTERING)
#error "ENABLE_ONEHOP_FILTERING requires ENABLE_SPOKE_FILTERING because OneHop consumes spoke-stage records."
#endif

// --- cde_edge_ie lower-bound options ---
#define CDE_LB_LIGHTWEIGHT_SPOKE

// --- other options ---
// #define NDEBUG

// --- cde_edge_ie anchor-support options ---
// Recompute anchor-support from the current search state when it is scored.
// This avoids expensive support maintenance in updateFrontier on high-result cases.
#define CDE_EDGE_IE_RECOMPUTE_ANCHOR_SUPPORT
// #define ENABLE_EXCLUDED_EDGE_SUPPORT

#endif //_CONFIG_H_
