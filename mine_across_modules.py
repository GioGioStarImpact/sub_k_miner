# -*- coding: utf-8 -*-
import sys, os, json, pandas as pd
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent))
import verilog_multi_to_aig as mm
from aig_subgraph_miner import summarize_patterns

def patterns_to_rows(patterns):
    rows = []
    items = sorted(patterns.items(), key=lambda kv: (-kv[1].count, kv[1].k, kv[1].m))
    for rank, (key, pat) in enumerate(items, 1):
        k, m, outs = key
        width = 1<<k
        def bits(tt, k): return [ (tt>>i)&1 for i in range(1<<k) ]
        cols = [bits(tt, k) for tt in outs]
        table_lines = []
        for i in range(width):
            in_bits = ' '.join(str((i>>b)&1) for b in range(k)) if k>0 else ''
            out_bits = ' '.join(str(cols[j][i]) for j in range(m)) if m>0 else ''
            table_lines.append(f"{i:3d} | {in_bits} || {out_bits}")
        rows.append({
            "rank": rank,
            "count": pat.count,
            "k": k,
            "m": m,
            "canon_outs_hash": f"{(hash(outs)&0xffffffff):08x}",
            "truth_table_2d": "\n".join(table_lines),
            "example": pat.example or {}
        })
    return rows

def patterns_to_matrix_rows(patterns):
    items = sorted(patterns.items(), key=lambda kv: (-kv[1].count, kv[1].k, kv[1].m))
    if not items:
        return [], 0, 0
    max_k = max(k for (k, m, outs) in (key for key, _ in items))
    max_m = max(m for (k, m, outs) in (key for key, _ in items))
    rows = []
    for rank, (key, pat) in enumerate(items, 1):
        k, m, outs = key
        width = 1 << k
        def bits(tt, k): return [ (tt>>i)&1 for i in range(1<<k) ]
        out_cols = [bits(tt, k) for tt in outs]
        for idx in range(width):
            row = {"rank": rank, "count": pat.count, "k": k, "m": m, "canon_outs_hash": f"{(hash(outs)&0xffffffff):08x}", "idx": idx}
            for b in range(max_k):
                row[f"x{b}"] = ((idx >> b) & 1) if b < k else ""
            for j in range(max_m):
                row[f"y{j}"] = out_cols[j][idx] if j < m else ""
            rows.append(row)
    return rows, max_k, max_m

def main():
    import argparse
    ap = argparse.ArgumentParser()
    ap.add_argument("verilog")
    ap.add_argument("--dmax", type=int, default=5)
    ap.add_argument("--kmax", type=int, default=8)
    ap.add_argument("--mmax", type=int, default=1)
    ap.add_argument("--min_support", type=int, default=1)
    ap.add_argument("--out_prefix", default="multi_module_report")
    ap.add_argument("--matrix_csv", action="store_true", help="Additionally write a matrix-friendly CSV with columns idx, x0.., y0..")
    args = ap.parse_args()

    vtext = Path(args.verilog).read_text(encoding="utf-8", errors="ignore")
    ir_per_module = mm.parse_verilog_to_ir_per_module(vtext)
    global_pats, per_mod = mm.mine_across_modules(ir_per_module, dmax=args.dmax, kmax=args.kmax, mmax=args.mmax, min_support=args.min_support)

    print("=== Global Summary ===")
    print(summarize_patterns(global_pats, top=50))
    for mod, pats in per_mod.items():
        print(f"\\n--- Module {mod} ---")
        print(summarize_patterns(pats, top=20))

    rows = patterns_to_rows(global_pats)
    import pandas as pd
    df = pd.DataFrame(rows)
    csv_path = f"/mnt/data/{args.out_prefix}_global.csv"
    json_path = f"/mnt/data/{args.out_prefix}_global.json"
    df.to_csv(csv_path, index=False)
    Path(json_path).write_text(json.dumps(rows, indent=2), encoding="utf-8")
    print(f"\\nSaved global CSV: {csv_path}")
    print(f"Saved global JSON: {json_path}")

    if args.matrix_csv:
        matrix_rows, max_k, max_m = patterns_to_matrix_rows(global_pats)
        mdf = pd.DataFrame(matrix_rows)
        matrix_csv_path = f"/mnt/data/{args.out_prefix}_global_matrix.csv"
        mdf.to_csv(matrix_csv_path, index=False)
        print(f"Saved matrix-friendly CSV: {matrix_csv_path}  (max_k={max_k}, max_m={max_m})")

if __name__ == "__main__":
    main()
