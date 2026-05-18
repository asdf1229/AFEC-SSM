#ifndef _CONFIG_H_
#define _CONFIG_H_

// --- cde_edge_ie filtering options ---
#define ENABLE_SPOKE_FILTERING
#define ENABLE_ONEHOP_FILTERING
#ifndef ONEHOP_FILTER_MISSING_GAP
#define ONEHOP_FILTER_MISSING_GAP 2
#endif
#define ENABLE_BRIDGE_FILTERING
// #define ENABLE_FRONTIER_ORDERING

// --- cde_edge_ie lower-bound options ---
#define CDE_LB_LIGHTWEIGHT_SPOKE

// --- other options ---
#define NDEBUG

#endif //_CONFIG_H_
