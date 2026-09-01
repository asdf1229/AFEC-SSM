#pragma once

// AFE-NoAF: keep AFE's filtering and root rule, but enumerate a
// branch-independent query-vertex order without anchor-frontier expansion.
#define AFEC_ENABLE_SPOKE_FILTERING 0
#define AFEC_ENABLE_BRIDGE_FILTERING 0
#define AFEC_ANCHOR_ORDER_FIXED 0
#define AFEC_ANCHOR_ORDER_DYNAMIC 1
#define AFEC_ANCHOR_ORDER_RANDOM 0
#define AFEC_BW_ALWAYS_BLACK 1
#define AFEC_BW_ALWAYS_WHITE 0
#define AFEC_BW_DYNAMIC 0
#define AFEC_ENABLE_COST_SPLIT 0
#define AFEC_ENABLE_ANCHOR_FRONTIER 0
