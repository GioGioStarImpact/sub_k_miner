# AIG Frequent Subgraph Miner (Vendor‑Neutral, Nangate-style cells)

This toolkit converts **gate-level Verilog** into **And‑Inverter Graphs (AIGs)** (one AIG per module), and mines
**frequent subgraphs** across the entire design. It includes:

- A lightweight multi‑module Verilog parser tuned for common Nangate45‑style cell names
  (`INV_X*`, `BUF_X*`, `NAND2_X*`, `NOR2_X*`, `DFF*`, `LOGIC0_X*`, `LOGIC1_X*`).
- An AIG core + frequent‑subgraph miner.
- A CLI runner that emits both human‑readable summaries (with **2D truth tables**) and
  matrix‑friendly CSV suitable for spreadsheets.

> Note: This project is **vendor‑neutral**. It does not contain or depend on any foundry‑specific content.

## Supported cell name patterns (heuristic)

- `INV_X*` → inverter (modeled as a buffer with inverted edge)
- `BUF_X*` → buffer (non‑inverting; if only `ZN` is present, treated as inverting buffer)
- `NAND2_X*` → 2‑input NAND (modeled as AND + inverted output)
- `NOR2_X*` → 2‑input NOR  (modeled as AND of inverted inputs)
- `DFF*` → flip‑flop cut‑points (treat `Q` as a **PI_from_seq** and `D` as a **PO_to_seq**; clock ignored)
- `LOGIC0_X*` / `LOGIC1_X*` → constant‑0 / constant‑1 sources (skipped by default to avoid polluting truth tables)

If your library uses different names, extend the simple name‑based rules in `verilog_multi_to_aig_nangate.py`.

## Installation

Python 3.9+ recommended.
```bash
pip install pandas
```

## Command‑line usage

```bash
python mine_across_modules_nangate.py /path/to/netlist.v \
  --dmax 4 --kmax 4 --mmax 1 --min_support 5 \
  --out_prefix run1 \
  --matrix_csv
```

**Outputs (in current directory):**
- `run1_ir_per_module.json` – per‑module IR dump (for inspection/debug)
- `run1_global.csv` – frequent patterns summary (each row is a canonical pattern; includes **2D truth table**)
- `run1_global.json` – same as JSON
- `run1_global_matrix.csv` – matrix‑friendly CSV (each input assignment is a row)

## Quick API usage

```python
from pathlib import Path
import verilog_multi_to_aig_nangate as mm
from aig_subgraph_miner import summarize_patterns

vtext = Path("netlist.v").read_text(encoding="utf-8", errors="ignore")
ir_per_module = mm.parse_verilog_to_ir_per_module(vtext)
global_pats, per_mod = mm.mine_across_modules(ir_per_module, dmax=4, kmax=4, mmax=1, min_support=5)

print("=== Global ===")
print(summarize_patterns(global_pats, top=20))
for mod, pats in per_mod.items():
    print(f"\n=== Module: {mod} ===")
    print(summarize_patterns(pats, top=10))
```

## Parameter guide

- `--dmax`: Maximum depth (levels) from cut leaves to root within a cone.
- `--kmax`: Maximum number of cut leaves (input variables) for truth tables.
- `--mmax`: Outputs per pattern (miner is currently single‑output `m=1`).
- `--min_support`: Global frequency threshold after aggregating across modules.

**Tips**
- Large designs: start conservative (e.g., `--kmax 4 --dmax 4 --min_support 5`), then relax gradually.
- Deeper mining: increase `--kmax` / `--dmax`, decrease `--min_support`.

## License

MIT – see `LICENSE`.
