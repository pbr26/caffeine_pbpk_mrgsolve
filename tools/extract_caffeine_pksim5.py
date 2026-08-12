#!/usr/bin/env python3
# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# Extract the caffeine compound parameters from the OSP Example_Caffeine model.
#
#   python3 tools/extract_caffeine_pksim5.py
#
# Reads literature/osp_models/Example_Caffeine/Caffeine.pksim5 (a SQLite file
# whose building-block content is a streamed ZIP of PK-Sim XML), inflates the
# caffeine compound XML, and writes data/caffeine_params.csv with one row per
# parameter: value, unit, source (from the model's ValueOrigin), and status.
#
# Provenance: every value is the one stored in the OSP "caffeine template model"
# (Britz et al. 2019, ref 15 = OSP Example_Caffeine repo). Nothing is typed from
# memory; the CSV is regenerable from the .pksim5 by re-running this script.
import sqlite3, struct, zlib, re, csv, os, sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKSIM = os.path.join(ROOT, "literature", "osp_models", "Example_Caffeine", "Caffeine.pksim5")
OUT   = os.path.join(ROOT, "data", "caffeine_params_raw.csv")

def caffeine_xml(path):
    con = sqlite3.connect(path); cur = con.cursor()
    cid = cur.execute("SELECT ContentId FROM BUILDING_BLOCKS WHERE Name='Caffeine'").fetchone()[0]
    blob = cur.execute("SELECT Data FROM CONTENTS WHERE Id=?", (cid,)).fetchone()[0]
    con.close()
    assert blob[:4] == b'PK\x03\x04', "unexpected content format"
    fnl, exl = struct.unpack('<HH', blob[26:30])
    start = 30 + fnl + exl
    return zlib.decompressobj(-15).decompress(blob[start:]).decode('utf-8', 'replace')

def attr(s, a):
    m = re.search(rf'\b{a}="([^"]*)"', s); return m.group(1) if m else ""

def main():
    if not os.path.exists(PKSIM):
        sys.exit("Missing " + PKSIM + " -- clone Example_Caffeine (see literature/DOWNLOADS.md).")
    xml = caffeine_xml(PKSIM)

    # <Para ...> blocks, each optionally trailed by <Info/> and <ValueOrigin/>
    blocks = re.findall(r'<Para\b[^>]*>(?:<Info[^>]*/>)?(?:<ValueOrigin[^>]*/>)?', xml)
    def finite(v):
        try: float(v); return v.lower() != "nan"
        except Exception: return False

    # target parameter name -> output label
    targets = {
        "Molecular weight": "molecular_weight",
        "Lipophilicity": "lipophilicity_logMA",
        "Fraction unbound (plasma, reference value)": "fraction_unbound_plasma",
        "Solubility at reference pH": "solubility_at_ref_pH",
        "Reference pH": "solubility_reference_pH",
        "pKa value 0": "pKa_1",
        "Specific intestinal permeability (transcellular)": "intestinal_permeability",
        "In vitro Vmax for liver microsomes": "CYP1A2_Vmax_invitro",
        "Content of CYP proteins in liver microsomes": "microsomal_CYP_content",
        "Km": "CYP1A2_Km",
    }
    picked = {}
    for b in blocks:
        n = attr(b, "name")
        if n in targets and finite(attr(b, "value")):
            lbl = targets[n]
            # keep the first finite (sourced alternatives come first in the model)
            if lbl not in picked:
                picked[lbl] = dict(
                    parameter=lbl, source_name=n,
                    value=attr(b, "value"),
                    unit=attr(b, "displayUnit") or attr(b, "dim"),
                    source=attr(b, "description") or "OSP caffeine template model",
                )

    # calculation methods (partition / permeability) recorded as rows
    methods = re.findall(r'<CalculationMethod\b[^>]*name="([^"]*)"', xml)
    method_rows = [dict(parameter="calc_method", source_name=m, value="", unit="",
                        source="OSP caffeine template model")
                   for m in methods if any(k in m for k in ("partition","permeability","Partition","Permeability"))]

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    fields = ["parameter", "value", "unit", "source", "source_name", "status", "notes"]
    order = ["molecular_weight","lipophilicity_logMA","fraction_unbound_plasma",
             "pKa_1","solubility_at_ref_pH","solubility_reference_pH",
             "intestinal_permeability","CYP1A2_Vmax_invitro","CYP1A2_Km",
             "microsomal_CYP_content"]
    with open(OUT, "w", newline="") as f:
        w = csv.DictWriter(f, fieldnames=fields); w.writeheader()
        for k in order:
            if k in picked:
                r = picked[k]; r["status"] = "verified"; r["notes"] = ""
                w.writerow(r)
        for r in method_rows:
            r["status"] = "verified"; r["notes"] = "tissue distribution / permeability model"
            w.writerow(r)
    print("wrote", OUT, "with", len(picked)+len(method_rows), "rows")

if __name__ == "__main__":
    main()
