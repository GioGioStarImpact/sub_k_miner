# -*- coding: utf-8 -*-
from __future__ import annotations
from typing import Dict, List, Tuple, Any
from pathlib import Path
import re, json

pin_re = re.compile(r"\.(\w+)\s*\(\s*([^()]+?)\s*\)")
inst_re = re.compile(r"(?ms)^\s*(\w+)\s+(\w+)\s*\((.*?)\)\s*;")

def _split_csv_names(s: str) -> List[str]:
    parts = [x.strip() for x in s.replace('\n',' ').split(',') if x.strip()]
    out = []
    for p in parts:
        p = re.sub(r"\[[^\]]+\]\s*", "", p).strip()
        if p: out.append(p)
    return out

def extract_modules(vtext: str) -> List[Tuple[str, str]]:
    mod_re = re.compile(r"(?ms)module\s+(\w+)\s*\(.*?\)\s*;(.+?)endmodule")
    mods = []
    for m in mod_re.finditer(vtext):
        name = m.group(1); body = m.group(2); mods.append((name, body))
    if not mods:
        m2 = re.search(r"module\s+(\w+)", vtext)
        name = m2.group(1) if m2 else "top"; mods = [(name, vtext)]
    return mods

def parse_module_to_ir(module_name: str, body_text: str) -> Dict[str, Any]:
    prim_in, prim_out = [], []
    for dirkw, bucket in (("input", prim_in), ("output", prim_out)):
        for m in re.finditer(rf"^\s*{dirkw}\b([^;]*);", body_text, flags=re.M):
            bucket.extend(_split_csv_names(m.group(1)))
    ands, bufs, pi_from_seq, po_to_seq, skipped = [], [], [], [], []
    for m in inst_re.finditer(body_text):
        cell, inst, pinblk = m.group(1), m.group(2), m.group(3)
        pins = {pm.group(1): pm.group(2).strip() for pm in re.finditer(pin_re, pinblk)}
        if cell == "INV":
            i = pins.get("I") or pins.get("A") or pins.get("IN")
            z = pins.get("ZN") or pins.get("Z") or pins.get("Y")
            if i and z: bufs.append({"out": z, "in": i, "inv": True})
            else: skipped.append({"type": cell, "inst": inst, "reason": "missing pins"})
        elif cell == "BUFF":
            i = pins.get("I") or pins.get("A") or pins.get("IN")
            z = pins.get("Z") or pins.get("ZN") or pins.get("Y")
            if i and z:
                inv = False
                if "ZN" in pins and "Z" not in pins: inv = True
                bufs.append({"out": z, "in": i, "inv": inv})
            else: skipped.append({"type": cell, "inst": inst, "reason": "missing pins"})
        elif cell == "ND2":
            a1 = pins.get("A1") or pins.get("A")
            a2 = pins.get("A2") or pins.get("B")
            z  = pins.get("ZN") or pins.get("Z")
            if a1 and a2 and z:
                if "ZN" in pins:
                    ands.append({"out": z, "in0": a1, "in0_inv": False, "in1": a2, "in1_inv": False, "out_inv": True})
                else:
                    ands.append({"out": z, "in0": a1, "in0_inv": False, "in1": a2, "in1_inv": False, "out_inv": False})
            else: skipped.append({"type": cell, "inst": inst, "reason": "missing pins"})
        elif cell == "NR2":
            a1 = pins.get("A1") or pins.get("A")
            a2 = pins.get("A2") or pins.get("B")
            z  = pins.get("ZN") or pins.get("Z")
            if a1 and a2 and z:
                if "ZN" in pins:
                    ands.append({"out": z, "in0": a1, "in0_inv": True, "in1": a2, "in1_inv": True, "out_inv": False})
                else:
                    ands.append({"out": z, "in0": a1, "in0_inv": True, "in1": a2, "in1_inv": True, "out_inv": True})
            else: skipped.append({"type": cell, "inst": inst, "reason": "missing pins"})
        elif cell in ("CKLNQ",):
            q = pins.get("Q"); e = pins.get("E") or pins.get("D")
            if q: pi_from_seq.append(q)
            if e: po_to_seq.append(e)
        elif cell in ("SDFQ","SDFRPQ","SDFSNQ"):
            q = pins.get("Q"); d = pins.get("D")
            if q: pi_from_seq.append(q)
            if d: po_to_seq.append(d)
        else:
            skipped.append({"type": cell, "inst": inst, "reason": "unsupported/blackbox"})
            continue
    return {"module": module_name, "primary_inputs": prim_in, "primary_outputs": prim_out, "ands": ands, "bufs": bufs, "pi_from_seq": sorted(set(pi_from_seq)), "po_to_seq": sorted(set(po_to_seq)), "skipped_cells": skipped}

def parse_verilog_to_ir_per_module(vtext: str) -> Dict[str, dict]:
    mods = extract_modules(vtext)
    return {name: parse_module_to_ir(name, body) for name, body in mods}

def build_aig_from_ir(ir: dict):
    from aig_subgraph_miner import AIG
    g = AIG(); net2ref = {}
    def ref_of(net: str):
        if net in net2ref: return net2ref[net]
        nid = g.add_pi(net); net2ref[net] = (nid, False); return net2ref[net]
    for n in ir.get("primary_inputs", []):
        if n not in net2ref: net2ref[n] = (g.add_pi(n), False)
    for n in ir.get("pi_from_seq", []):
        if n not in net2ref: net2ref[n] = (g.add_pi(n), False)
    for b in ir.get("bufs", []):
        s_id, s_inv = ref_of(b["in"]); net2ref[b["out"]] = (s_id, s_inv ^ bool(b.get("inv", False)))
    for a in ir.get("ands", []):
        in0_id, in0_inv = ref_of(a["in0"]); in0_inv ^= bool(a.get("in0_inv", False))
        in1_id, in1_inv = ref_of(a["in1"]); in1_inv ^= bool(a.get("in1_inv", False))
        nid = g.add_and(in0_id, in1_id, inv_a=in0_inv, inv_b=in1_inv, name=a.get("out"))
        net2ref[a["out"]] = (nid, bool(a.get("out_inv", False)))
    for name in (ir.get("primary_outputs", []) + ir.get("po_to_seq", [])):
        if name in net2ref:
            nid, inv = net2ref[name]; g.add_po(nid, inv=inv, name=name)
        else:
            nid = g.add_pi(name+"__dangling"); g.add_po(nid, inv=False, name=name)
    return g

def mine_across_modules(ir_per_module: Dict[str, dict], dmax=5, kmax=8, mmax=1, min_support=1):
    from aig_subgraph_miner import mine_frequent_subgraphs, Pattern
    global_patterns = {}
    per_module_patterns = {}
    for mod, ir in ir_per_module.items():
        g = build_aig_from_ir(ir)
        pats = mine_frequent_subgraphs(g, dmax=dmax, kmax=kmax, mmax=mmax, min_support=1)
        per_module_patterns[mod] = pats
        for key, p in pats.items():
            if key not in global_patterns:
                global_patterns[key] = Pattern(k=p.k, m=p.m, canonical_outputs=p.canonical_outputs, count=p.count, example=p.example)
            else:
                global_patterns[key].count += p.count
    global_patterns = {k:v for k,v in global_patterns.items() if v.count >= min_support}
    return global_patterns, per_module_patterns
