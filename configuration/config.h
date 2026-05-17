#ifndef _CONFIG_H_
#define _CONFIG_H_

// --- filtering options ---
// #define ENABLE_EDGE_LABEL
#define ENABLE_ADVANCED_FILTERING
#define ENABLE_SPOKE_FILTERING
#define ENABLE_ONEHOP_FILTERING
#ifndef ONEHOP_FILTER_MISSING_GAP
#define ONEHOP_FILTER_MISSING_GAP 2
#endif
#define ENABLE_BRIDGE_FILTERING
// #define ENABLE_FRONTIER_ORDERING

// --- lower-bound options ---
// #define LOWER_BOUND
#define CDE_LB_LIGHTWEIGHT_SPOKE
// #define CDE_LB_COMPONENT_MWPM
// #define CDE_LB_COMPONENT_MWPM_CACHE
#if defined(CDE_LB_LIGHTWEIGHT_SPOKE)
#define LOWER_BOUND
#endif
#ifndef LOWER_BOUND_MISSING_GAP
#define LOWER_BOUND_MISSING_GAP 0
#endif

// --- other options ---
// #define NDEBUG

#endif //_CONFIG_H_
