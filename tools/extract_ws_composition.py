#!/usr/bin/env python3
# Author: Pramod B R  |  caffeine_pbpk_mrgsolve
# Extract PK-Sim's Willmann-Schmitt (-WS) tissue composition (fractional volume
# of water / neutral lipid / neutral phospholipid per organ) from a PK-Sim .pkml
# simulation. Composition is drug-independent physiology, so ANY human .pkml
# works -- we use ospsuite's bundled Aciclovir.pkml (from the R package extdata).
#
# How the .pkml was obtained (fully reproducible in R):
#   pkmls <- list.files(system.file("extdata", package="ospsuite"),
#                       pattern="\\.pkml$", full.names=TRUE)
#   file.copy(pkmls[1], "tables/Aciclovir.pkml")   # the Aciclovir human sim
# Fat & Muscle are individual-specific (not inline in the .pkml); those two rows
# come from ospsuite::createIndividual(European ICRP 2002) instead.
#
#   python3 tools/extract_ws_composition.py tables/Aciclovir.pkml
import sys, xml.etree.ElementTree as ET
ORGANS={"Bone","Brain","Fat","Gonads","Heart","Kidney","Liver","Lung","Muscle",
        "Skin","Stomach","SmallIntestine","LargeIntestine","Pancreas","Spleen"}
WANT={"Vf (water)-WS","Vf (neutral lipid)-WS","Vf (neutral phospholipid, plasma)-WS"}
def main(path):
    root=ET.parse(path).getroot(); res={}
    def walk(el,organ):
        nm=el.get("name")
        if el.tag.endswith("Container") and nm in ORGANS: organ=nm
        if el.tag.endswith("Parameter") and nm in WANT and el.get("value") is not None:
            res.setdefault(organ,{})[nm]=el.get("value")
        for ch in el: walk(ch,organ)
    walk(root,None)
    print("organ,Vw,Vnl,Vph")
    for org in sorted(res):
        w=res[org]
        print(f"{org},{w.get('Vf (water)-WS','')},{w.get('Vf (neutral lipid)-WS','')},"
              f"{w.get('Vf (neutral phospholipid, plasma)-WS','')}")
if __name__=="__main__":
    main(sys.argv[1] if len(sys.argv)>1 else "tables/Aciclovir.pkml")
