#pragma once

#include <cstddef>

#include <sdsl/rank_support_v.hpp>
#include <sdsl/select_support_mcl.hpp>
#include <sdsl/wm_int.hpp>

namespace lz77tax {

struct SharedWm : public sdsl::wm_int<> {
    size_t zero_cnt_at(size_t k)   const { return m_zero_cnt[k]; }
    size_t rank_level_at(size_t k) const { return m_rank_level[k]; }
    size_t tree_rank1(size_t i)    const { return m_tree_rank(i); }
    size_t tree_sel0(size_t i)     const { return m_tree_select0(i); }
    size_t tree_sel1(size_t i)     const { return m_tree_select1(i); }
};

}  // namespace lz77tax
