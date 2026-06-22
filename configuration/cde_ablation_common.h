#pragma once

// Shared settings for every CDE-Edge-IE ablation executable.
// Define a value in a variant header before including this file to override it.

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

#ifndef CDE_EDGE_IE_FIXED_ORDER
#define CDE_EDGE_IE_FIXED_ORDER 0
#endif

#define NDEBUG
