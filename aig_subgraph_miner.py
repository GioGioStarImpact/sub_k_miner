# -*- coding: utf-8 -*-
from dataclasses import dataclass
from typing import *
from collections import defaultdict, deque
import itertools, random

@dataclass(frozen=True)
class Edge:
    src: int
    inv: bool = False

@dataclass
class Node:
    nid: int
    kind: str
    fin0: Optional[Edge] = None
    fin1: Optional[Edge] = None
    name: Optional[str] = None

class AIG:
    def __init__(self):
        self.nodes: Dict[int, Node] = {}
        self._fanouts: Dict[int, list] = defaultdict(list)
        self._pis: list = []
        self._pos: list = []
        self._next_id: int = 1

    def add_pi(self, name: Optional[str] = None) -> int:
        nid = self._next_id; self._next_id += 1
        self.nodes[nid] = Node(nid, 'PI', name=name)
        self._pis.append(nid); return nid

    def add_and(self, a: int, b: int, inv_a: bool=False, inv_b: bool=False, name: Optional[str]=None) -> int:
        nid = self._next_id; self._next_id += 1
        n = Node(nid, 'AND', Edge(a, inv_a), Edge(b, inv_b), name)
        self.nodes[nid] = n
        self._fanouts[a].append(Edge(nid)); self._fanouts[b].append(Edge(nid))
        return nid

    def add_po(self, src: int, inv: bool=False, name: Optional[str]=None) -> int:
        nid = self._next_id; self._next_id += 1
        n = Node(nid, 'PO', Edge(src, inv), name=name)
        self.nodes[nid] = n
        self._fanouts[src].append(Edge(nid)); self._pos.append(nid)
        return nid

    def fanouts(self, nid: int): return self._fanouts.get(nid, [])

    def topological_order(self) -> list:
        indeg = defaultdict(int)
        for n in self.nodes.values():
            if n.kind == 'AND':
                indeg[n.nid] += 2; indeg[n.fin0.src] += 0; indeg[n.fin1.src] += 0
            elif n.kind == 'PO':
                indeg[n.nid] += 1; indeg[n.fin0.src] += 0
            else:
                indeg[n.nid] += 0
        q = deque([nid for nid, d in indeg.items() if d == 0]); order = []
        while q:
            u = q.popleft(); order.append(u)
            for e in self.fanouts(u):
                indeg[e.src] -= 1
                if indeg[e.src] == 0: q.append(e.src)
        seen = set(order)
        for n in self.nodes:
            if n not in seen: order.append(n)
        return order

    def compute_levels(self) -> Dict[int,int]:
        level: Dict[int,int] = {}
        for nid in self.topological_order():
            n = self.nodes[nid]
            if n.kind == 'PI': level[nid] = 0
            elif n.kind == 'AND': level[nid] = 1 + max(level[n.fin0.src], level[n.fin1.src])
            elif n.kind == 'PO': level[nid] = level[n.fin0.src]
        return level

    def cone_nodes(self, roots: Iterable[int], stop_at: set) -> list:
        roots = list(roots); visited = set(); stack = list(roots)
        while stack:
            u = stack.pop()
            if u in visited: continue
            visited.add(u); n = self.nodes[u]
            if n.kind == 'AND':
                for e in (n.fin0, n.fin1):
                    if e.src not in visited and e.src not in stop_at: stack.append(e.src)
            elif n.kind == 'PO':
                if n.fin0.src not in visited and n.fin0.src not in stop_at: stack.append(n.fin0.src)
        full = self.topological_order()
        return [nid for nid in full if nid in visited]

def enumerate_k_feasible_cuts(aig, root, kmax, dmax, levels, cap_per_node=600):
    memo = {}
    def _cuts(nid):
        if nid in memo: return memo[nid]
        n = aig.nodes[nid]; out = set()
        if n.kind == 'PI':
            out.add(frozenset([nid]))
        elif n.kind == 'AND':
            c0 = _cuts(n.fin0.src); c1 = _cuts(n.fin1.src)
            for a in c0:
                for b in c1:
                    m = a | b
                    if len(m) <= kmax and all((levels[nid]-levels[l]) < dmax for l in m):
                        out.add(m)
                if len(out) > cap_per_node: break
        elif n.kind == 'PO':
            out = _cuts(n.fin0.src)
        memo[nid] = out; return out
    cuts = [tuple(sorted(fs)) for fs in _cuts(root)]
    cuts.sort(key=lambda t: (len(t), t)); return cuts

def var_patterns(k):
    n = 1 << k; pats = []
    for i in range(k):
        block = 1 << i; pattern = 0; bit = 0
        while bit < n:
            bit += block
            for _ in range(block):
                if bit >= n: break
                pattern |= (1 << bit); bit += 1
        pats.append(pattern)
    return pats

def invert_bits(x, width): return (~x) & ((1 << width) - 1)

def permute_truth_table(tt, k, perm, phases):
    if k == 0: return tt
    width = 1 << k; out = 0
    for new_idx in range(width):
        orig_idx = 0
        for new_pos, orig_var in enumerate(perm):
            bit = (new_idx >> new_pos) & 1
            if (phases >> orig_var) & 1: bit ^= 1
            if bit: orig_idx |= (1 << orig_var)
        if (tt >> orig_idx) & 1: out |= (1 << new_idx)
    return out

def tt_lex_maximize_over_output_phase(tts, k):
    width = 1 << k; outs = []; phases = []
    for f in tts:
        neg = invert_bits(f, width)
        if f >= neg: outs.append(f); phases.append(False)
        else: outs.append(neg); phases.append(True)
    return tuple(outs), tuple(phases)

def literal_signature(tts, k, var):
    n = 1 << k; ones0 = ones1 = 0; mask = 1 << var
    for idx in range(n):
        if (idx & mask) == 0:
            for f in tts: ones0 += (f >> idx) & 1
        else:
            for f in tts: ones1 += (f >> idx) & 1
    return (ones0, ones1)

def group_variables_by_signature(tts, k):
    from collections import defaultdict
    sig2vars = defaultdict(list)
    for v in range(k):
        sig2vars[literal_signature(tts, k, v)].append(v)
    return [sorted(vs) for vs in sig2vars.values()]

def npn_canon_multi_output(tts, k, perm_cap=120000):
    if k == 0:
        outs, outph = tt_lex_maximize_over_output_phase(tts, k)
        outs = tuple(sorted(outs, reverse=True))
        return outs, tuple(), 0, outph
    groups = group_variables_by_signature(tts, k)
    group_perms = [list(itertools.permutations(g)) for g in groups]
    best_outs_sorted = None; best_perm = None; best_phases = None; best_outph = None
    for local in itertools.product(*group_perms):
        perm_list = [None]*k; pos = 0
        for seq in local:
            for v in seq:
                perm_list[pos] = v; pos += 1
        perm = tuple(perm_list)
        phases = 0  # simplification here
        outs_t = tuple(permute_truth_table(f, k, perm, phases) for f in tts)
        outs_sorted, outph = tt_lex_maximize_over_output_phase(outs_t, k)
        outs_sorted = tuple(sorted(outs_sorted, reverse=True))
        if (best_outs_sorted is None) or (outs_sorted > best_outs_sorted):
            best_outs_sorted, best_perm, best_phases, best_outph = outs_sorted, perm, phases, outph
    return best_outs_sorted, best_perm, best_phases, best_outph

def random_sim_signature(aig, roots, leaves, num_vecs=64, seed=1):
    rng = random.Random(seed + sum(roots) + sum(leaves)*1315423911)
    cone = aig.cone_nodes(roots, stop_at=set(leaves))
    sig = 0
    for _ in range(num_vecs):
        val_leaf = {leaf: rng.randint(0,1) for leaf in leaves}; val = {}
        for nid in cone:
            n = aig.nodes[nid]
            if n.kind == 'PI':
                val[nid] = val_leaf.get(nid, 0)
            elif n.kind == 'AND':
                a = val.get(n.fin0.src, val_leaf.get(n.fin0.src)); b = val.get(n.fin1.src, val_leaf.get(n.fin1.src))
                a ^= int(n.fin0.inv); b ^= int(n.fin1.inv); val[nid] = a & b
            elif n.kind == 'PO':
                x = val.get(n.fin0.src, val_leaf.get(n.fin0.src)); val[nid] = x ^ int(n.fin0.inv)
        outbits = 0
        for i, r in enumerate(roots): outbits |= (val[r] & 1) << i
        sig = ((sig << len(roots)) | outbits) & ((1 << 256) - 1)
    return sig

def eval_truth_tables(aig, roots, leaves):
    k = len(leaves); nbits = 1 << k; pats = var_patterns(k)
    leaf_tt = {leaf: pats[i] for i, leaf in enumerate(leaves)}
    cone = aig.cone_nodes(roots, stop_at=set(leaves)); val = {}
    for nid in cone:
        n = aig.nodes[nid]
        if n.kind == 'PI':
            if n.nid in leaf_tt: val[nid] = leaf_tt[nid]
        elif n.kind == 'AND':
            a = val.get(n.fin0.src, leaf_tt.get(n.fin0.src)); b = val.get(n.fin1.src, leaf_tt.get(n.fin1.src))
            if n.fin0.inv: a = invert_bits(a, nbits)
            if n.fin1.inv: b = invert_bits(b, nbits)
            val[nid] = a & b
        elif n.kind == 'PO':
            x = val.get(n.fin0.src, leaf_tt.get(n.fin0.src))
            if n.fin0.inv: x = invert_bits(x, nbits)
            val[nid] = x
    return tuple(val[r] for r in roots)

@dataclass
class Pattern:
    k: int; m: int; canonical_outputs: Tuple[int,...]; count: int = 0; example: Optional[Dict] = None

def mine_frequent_subgraphs(aig, dmax=5, kmax=8, mmax=1, min_support=1, cap_per_node=600, use_random_bucket=True):
    levels = aig.compute_levels(); order = aig.topological_order()
    roots = [nid for nid in order if aig.nodes[nid].kind in ('AND','PO')]
    patterns: Dict[Tuple[int,int,Tuple[int,...]], Pattern] = {}
    canon_cache = {}
    for r in roots:
        cuts = enumerate_k_feasible_cuts(aig, r, kmax, dmax, levels, cap_per_node=cap_per_node)
        for leaves in cuts:
            k = len(leaves)
            if k == 0 or k > kmax: continue
            sig = random_sim_signature(aig, (r,), leaves, 64, 7) if use_random_bucket else 0
            cache_key = (k, 1, hash(sig) & 0xffffffff)
            outs = eval_truth_tables(aig, (r,), leaves)
            if use_random_bucket and cache_key in canon_cache:
                canon_outs, perm, inph, outph = canon_cache[cache_key]
            else:
                canon_outs, perm, inph, outph = npn_canon_multi_output(outs, k)
                if use_random_bucket: canon_cache[cache_key] = (canon_outs, perm, inph, outph)
            key = (k, 1, canon_outs)
            p = patterns.get(key)
            if p is None:
                patterns[key] = Pattern(k=k, m=1, canonical_outputs=canon_outs, count=1,
                    example={"root": r, "leaves": leaves, "perm_newpos_to_origvar": perm, "input_phase_mask_on_origvars": inph, "output_phases": outph})
            else:
                p.count += 1
    return {k:v for k,v in patterns.items() if v.count >= min_support}

def summarize_patterns(patterns, top: int = 20) -> str:
    def _tt_bits(tt: int, k: int) -> list:
        width = 1 << k
        return [(tt >> i) & 1 for i in range(width)]
    items = sorted(patterns.items(), key=lambda kv: (-kv[1].count, kv[1].k, kv[1].m))
    lines: list = []
    for idx, (key, pat) in enumerate(items[:top], 1):
        k, m, outs = key
        lines.append(f"[{idx:02d}] count={pat.count}  k={k}  m={m}  canon_outs_hash={(hash(outs) & 0xffffffff):08x}")
        width = 1 << k if k >= 0 else 0
        in_header = " ".join([f"x{j}" for j in range(k)]) if k > 0 else ""
        out_header = " ".join([f"y{j}" for j in range(m)]) if m > 0 else ""
        lines.append(f"    Truth table (rows: 0..{max(0, width-1)})")
        lines.append(f"    idx | {in_header:^{max(1, 2*k-1)}} || {out_header}")
        lines.append(f"    {'-'*3}-+-{'-'*max(1, 2*k-1)}-++-{'-'*max(3, 2*m-1)}")
        cols = [_tt_bits(tt, k) for tt in outs]
        for i in range(width):
            bits_in = [(i >> b) & 1 for b in range(k)]
            in_str = " ".join(str(b) for b in bits_in) if k > 0 else ""
            out_str = " ".join(str(cols[j][i]) for j in range(m)) if m > 0 else ""
            lines.append(f"    {i:3d} | {in_str:^{max(1, 2*k-1)}} || {out_str}")
    return "\n".join(lines)
