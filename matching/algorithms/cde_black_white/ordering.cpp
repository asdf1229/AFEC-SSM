#include "matching/algorithms/cde_black_white/cde_black_white.h"

namespace cde_black_white {

ui MatchingSolver::chooseRoot()
{
    ui root = 0;
    for (ui u = 1; u < qn; ++u) {
        size_t cand_u = candidates[u].size();
        size_t cand_root = candidates[root].size();

        ui deg_u = q_degree[u];
        ui deg_root = q_degree[root];

        if (cand_u * deg_root < cand_root * deg_u) root = u;
    }
    return root;
}

void MatchingSolver::initColors()
{
    static_root = chooseRoot();
    static_color.assign(qn, COLOR_WHITE);
    static_color[static_root] = COLOR_BLACK;
}

ui MatchingSolver::chooseMatWhite(const SearchState &state) const
{
    // 选择候选数最少的 white 点作为优先具体化对象。
    ui chosen = qn;
    ui best_count = std::numeric_limits<ui>::max();
    for (ui u = 0; u < qn; ++u) {
        if (!isWhite(state, u)) {
            continue;
        }
        ui count = state.white[u].feasible_count;
        if (count < best_count || (count == best_count && u < chosen)) {
            chosen = u;
            best_count = count;
        }
    }
    return chosen;
}

} // namespace cde_black_white
