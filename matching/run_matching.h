#ifndef SSM_MATCHING_RUN_MATCHING_H_
#define SSM_MATCHING_RUN_MATCHING_H_

#include "graph/graph.h"

#include <limits>

namespace ssm {

    using MatchResults = std::vector<std::vector<std::pair<ui, ui> > >;
    using AlgorithmEntry = void (*)(const Graph *, const Graph *, MatchResults &, ui);

    struct AlgorithmDefinition {
        std::string key;
        std::string display_name;
        AlgorithmEntry entry;
    };

    void clear_reported_result_count();
    void set_reported_result_count(size_t count);
    size_t get_reported_result_count(const MatchResults &results);

    namespace detail {
        extern size_t next_output_report_count;
        void emit_output_checkpoint(size_t count, unsigned long long recursion_calls);
    }

    // Keep the per-result hot path to one predictable comparison.  The slow
    // path only runs at the configured output checkpoints.
    inline void report_result_progress(
        size_t count,
        unsigned long long recursion_calls = std::numeric_limits<unsigned long long>::max())
    {
        if (count >= detail::next_output_report_count) {
            detail::emit_output_checkpoint(count, recursion_calls);
        }
    }

    size_t peak_resident_memory_kb();

    struct AlgorithmPhaseTimes {
        long long preprocessing_us = 0;
        long long search_us = 0;
        bool available = false;
    };

    void clear_reported_phase_times();
    void report_preprocessing_complete(long long preprocessing_us);
    void set_reported_phase_times(long long preprocessing_us, long long search_us);
    AlgorithmPhaseTimes get_reported_phase_times();

    struct AlgorithmKeyMetrics {
        size_t filter_candidates = 0;
        unsigned long long recursion_calls = 0;
        bool filter_candidates_available = false;
        bool recursion_calls_available = false;
    };

    void clear_reported_key_metrics();
    void set_reported_filter_candidates(size_t count);
    void set_reported_recursion_calls(unsigned long long count);
    AlgorithmKeyMetrics get_reported_key_metrics();

    const AlgorithmDefinition &create_algorithm_definition();
    int run_algorithm_main(int argc, char *argv[], const AlgorithmDefinition &algorithm);

} // namespace ssm

#endif
