#!/usr/bin/env python3
# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# Extract the observed caffeine plasma profiles from the OSP
# Fluvoxamine-Caffeine-DDI PK-Sim project into data/digitized/.
#
#   python3 tools/extract_observed_pksim5.py
#
# The observed data are stored in the .pksim5 (SQLite) as building-block
# CONTENTS: a streamed ZIP of PK-Sim XML, whose <Values> arrays are base64
# .NET-serialized float32 vectors in PK-Sim BASE units (time = min,
# mass concentration = kg/L). This script inflates, decodes, converts to
# display units (h, ug/mL, ng/mL) and time-after-dose, and writes one CSV per
# study. Regenerable from the repo; nothing typed by hand.
#
# Source: Britz et al. 2019 (OSP Fluvoxamine-Caffeine-DDI). Underlying clinical
# studies: Culm-Merdek et al. 2005 and Jeppesen et al. 1996.
import sqlite3, struct, zlib, re, base64, csv, os, sys

ROOT  = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PKSIM = os.path.join(ROOT, "literature", "osp_models",
                     "Fluvoxamine-Caffeine-DDI", "Fluvoxamine-Caffeine-DDI.pksim5")
OUT   = os.path.join(ROOT, "data", "digitized")

# caffeine "alone" (no fluvoxamine) arms -> the profiles to validate the base model
STUDIES = {
    "CulmMerdek2005_250mgSD_NA_NA_Caffeine_Alone": ("culm_merdek_2005_250mg.csv",
                                                    "Culm-Merdek et al. 2005, caffeine 250 mg single oral dose"),
    "Jeppesen1996_200mgSD_EM_NO_Caffeine_Alone":   ("jeppesen_1996_200mg.csv",
                                                    "Jeppesen et al. 1996, caffeine 200 mg single oral dose (EM)"),
}

def inflate(con, cid):
    b = con.execute("SELECT Data FROM CONTENTS WHERE Id=?", (cid,)).fetchone()[0]
    fnl, exl = struct.unpack('<HH', b[26:30]); start = 30 + fnl + exl
    return zlib.decompressobj(-15).decompress(b[start:]).decode('utf-8', 'replace')

def floats(b64):
    raw = base64.b64decode(b64)
    n = struct.unpack('<i', raw[22:26])[0]          # ArraySinglePrimitive length
    return list(struct.unpack('<%df' % n, raw[27:27 + 4 * n]))

def main():
    if not os.path.exists(PKSIM):
        sys.exit("Missing " + PKSIM + " (clone Fluvoxamine-Caffeine-DDI; see literature/DOWNLOADS.md).")
    os.makedirs(OUT, exist_ok=True)
    con = sqlite3.connect(PKSIM)
    for name, (fname, desc) in STUDIES.items():
        cid = con.execute("SELECT ContentId FROM DATA_REPOSITORIES WHERE Name=?", (name,)).fetchone()[0]
        xml = inflate(con, cid)
        blocks = re.findall(r'<Values>([^<]+)</Values>', xml)   # order: Time, Conc, SD
        t, c, sd = floats(blocks[0]), floats(blocks[1]), floats(blocks[2])
        t0 = t[0]
        path = os.path.join(OUT, fname)
        with open(path, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["time_h", "conc_ug_mL", "sd_ng_mL", "study", "source"])
            for ti, ci, si in zip(t, c, sd):
                w.writerow([round((ti - t0) / 60, 4),        # min -> h, time after dose
                            round(ci * 1e6, 5),              # kg/L -> ug/mL
                            round(si * 1e9, 3),              # kg/L -> ng/mL
                            desc, "OSP Fluvoxamine-Caffeine-DDI (Britz 2019)"])
        print("wrote", path, f"({len(t)} points, Cmax {max(c)*1e6:.2f} ug/mL)")
    con.close()

if __name__ == "__main__":
    main()
