#include "matching/algorithms/cde_black_white/cde_black_white.h"

namespace cde_black_white {

bool MatchingSolver::tryBindRoot(SearchState &state, ui root, ui v) const
{
    // 尝试把根查询点绑定到数据点 v，作为搜索初始 black 映射。
    if (root >= qn || v >= gn || !candidates[root].contains(v)) {
        return false;
    }
    state.color[root] = COLOR_BLACK;
    state.mapped_q[root] = (int)v;
    state.used_data_vertices.push_back(v);
    state.used_data_flag[v] = 1;
    state.part_M.push_back({ root, v });
    state.selected_count = 1;
    refreshFrontierEdgesIncidentTo(state, root);
    return true;
}

} // namespace cde_black_white
