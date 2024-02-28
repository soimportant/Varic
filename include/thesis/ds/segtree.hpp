#pragma once

#include <algorithm>
#include <cassert>
#include <cstdint>
#include <functional>
#include <vector>

/**
 * We use segment tree to define the minimum coverage that a set of intervals
 * can cover when a read has enough coverage, the we can say that the read is
 * valid for assembly
 */

class segtree {
 public:
  class node {
   public:
    /**
     * Initialize the data member of leaves
     */
    node() {
      min_coverage = 0;
      max_coverage = 0;
      lazy_add = 0;
    }

    /**
     * determine how to modify the property of a segment
     * remember to set the lazy flag
     */
    void apply(int val) {
      min_coverage += val;
      max_coverage += val;
      lazy_add += val;
    }
    std::int16_t min_coverage;
    std::int16_t max_coverage;
    std::int16_t lazy_add = 0;
  };

  segtree(int _n) : n(_n) {
    assert(n > 0);
    int sz = 1;
    while (sz < n) {
      sz <<= 1;
    }
    tree.resize(sz << 1);
    build(1, 0, n - 1);
  }

  /* move constructor */
  segtree(segtree&& other) {
    n = other.n;
    tree = std::move(other.tree);
  }

  /* move assignment */
  segtree& operator=(segtree&& other) {
    n = other.n;
    tree = std::move(other.tree);
    return *this;
  }

  /**
   * get the merge result from two segments
   */
  node merge(const node& a, const node& b) {
    node res{};
    res.min_coverage = std::min(a.min_coverage, b.min_coverage);
    res.max_coverage = std::max(a.max_coverage, b.max_coverage);
    return res;
  }

  /**
   * push property of the node to its child node
   */
  void push(int idx) {
    int idxL = idx * 2, idxR = idx * 2 + 1;

    /**
     * if the lazy flag is set, push the property to its child
     */
    if (tree[idx].lazy_add != 0) {
      tree[idxL].apply(tree[idx].lazy_add);
      tree[idxR].apply(tree[idx].lazy_add);
      tree[idx].lazy_add = 0;
    }
  }

  // template<class T>
  // segtree(const std::vector<T>& v) {
  //   n = (int) v.size();
  //   assert(n > 0);
  //   int sz = 1;
  //   while (sz < n) {
  //     sz <<= 1;
  //   }
  //   tree.resize(sz << 1);
  //   build(v, 1, 0, n - 1);
  // }

  // template<class T>
  // void build(const std::vector<T>& v) {
  //   assert((int) v.size() == n);
  //   build(v, 1, 0, n - 1);
  // }

  template <class... T>
  void modify(int l, int r, const T&... val) {
    assert(0 <= l && l <= r && r < n);
    modify(l, r, 1, 0, n - 1, val...);
  }
  node get(int l, int r) {
    assert(0 <= l && l <= r && r < n);
    return get(l, r, 1, 0, n - 1);
  }
  node get(int p) {
    assert(0 <= p && p < n);
    return get(p, p, 1, 0, n - 1);
  }

  // int find_last(int l, int r, const function<bool(const node&)>& f) {
  //   return find_last(l, r, 1, 0, n - 1, f);
  // }

  // int find_first(int l, int r, const function<bool(const node&)>& f) {
  //   return find_first(l, r, 1, 0, n - 1, f);
  // }

 private:
  void build(int idx, int l, int r) {
    if (l == r) {
      return;
    }
    int mid = (l + r) >> 1, idxL = idx * 2, idxR = idx * 2 + 1;
    build(idxL, l, mid);
    build(idxR, mid + 1, r);
    tree[idx] = merge(tree[idxL], tree[idxR]);
  }

  // template<class T>
  // void build(const std::vector<T>& v, int idx, int l, int r) {
  //   if (l == r) {
  //     tree[idx].apply(v[l]);
  //     return;
  //   }
  //   int mid = (l + r) >> 1, idxL = idx * 2, idxR = idx * 2 + 1;
  //   build(v, idxL, l, mid);
  //   build(v, idxR, mid + 1, r);
  //   tree[idx] = merge(tree[idxL], tree[idxR]);
  // }

  template <class... T>
  void modify(const int& ql, const int& qr, int idx, int l, int r, T... val) {
    if (ql <= l && r <= qr) {
      tree[idx].apply(val...);
      return;
    }
    push(idx);
    int mid = (l + r) >> 1, idxL = idx * 2, idxR = idx * 2 + 1;
    if (ql <= mid) {
      modify(ql, qr, idxL, l, mid, val...);
    }
    if (qr > mid) {
      modify(ql, qr, idxR, mid + 1, r, val...);
    }
    tree[idx] = merge(tree[idxL], tree[idxR]);
  }

  node get(int ql, int qr, int idx, int l, int r) {
    if (ql <= l && r <= qr) {
      return tree[idx];
    }
    push(idx);
    int mid = (l + r) >> 1, idxL = idx * 2, idxR = idx * 2 + 1;
    node res{};
    if (qr <= mid) {
      res = get(ql, qr, idxL, l, mid);
    } else if (ql > mid) {
      res = get(ql, qr, idxR, mid + 1, r);
    } else {
      res = merge(get(ql, qr, idxL, l, mid), get(ql, qr, idxR, mid + 1, r));
    }
    return res;
  }

  // int find_last(int ql, int qr, int idx, int l, int r, const
  // function<bool(const node&)>& f) {
  //   int mid = (l + r) >> 1, idxL = idx * 2, idxR = idx * 2 + 1;
  //   if (l > qr || r < ql) {
  //     return -1;
  //   }
  //   if (!f(tree[idx])) {
  //     return -1;
  //   } else {
  //     if (l == r) {
  //       return l;
  //     }
  //   }
  //   push(idx);
  //   int res = find_last(ql, qr, idxR, mid + 1, r, f);
  //   if (res == -1) {
  //     res = find_last(ql, qr, idxL, l, mid, f);
  //   }
  //   return res;
  // }

  // int find_first(int ql, int qr, int idx, int l, int r, const
  // function<bool(const node&)>& f) {
  //   int mid = (l + r) >> 1, idxL = idx * 2, idxR = idx * 2 + 1;
  //   if (l > qr || r < ql) {
  //     return -1;
  //   }
  //   if (!f(tree[idx])) {
  //     return -1;
  //   } else {
  //     if (l == r) {
  //       return l;
  //     }
  //   }
  //   push(idx);
  //   int res = find_first(ql, qr, idxL, l, mid, f);
  //   if (res == -1) {
  //     res = find_first(ql, qr, idxR, mid + 1, r, f);
  //   }
  //   return res;
  // }

  int n;
  std::vector<node> tree;
};

using segnode = segtree::node;
