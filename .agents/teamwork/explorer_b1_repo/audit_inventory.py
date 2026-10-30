import hashlib
import os
from pathlib import Path
import pandas as pd

root = Path(r"C:\Users\hp\.gemini\antigravity\scratch\telomerase-research")

def sha256_file(p: Path) -> str:
    h = hashlib.sha256()
    with open(p, "rb") as f:
        while chunk := f.read(65536):
            h.update(chunk)
    return h.hexdigest().lower()

print("--- 1. Raw Data Checksums & Manifest Comparison ---")
manifest_path = root / "data_manifest.csv"
df_man = pd.read_csv(manifest_path)
for _, r in df_man.iterrows():
    pid = r["pdb_id"]
    exp_h = r["sha256"]
    cif_p = root / "data" / "raw" / f"{pid}.cif"
    if cif_p.exists():
        act_h = sha256_file(cif_p)
        match = (exp_h == act_h)
        sz = cif_p.stat().st_size
        print(f"PDB {pid}: size={sz} bytes, match={match}, hash={act_h[:12]}...")
    else:
        print(f"PDB {pid}: MISSING!")

print("\n--- 2. Results Tables Inventory ---")
tables_dir = root / "results" / "tables"
for f in sorted(tables_dir.glob("*.csv")):
    df = pd.read_csv(f)
    print(f"{f.name:<35} | rows={len(df):<4} | cols={len(df.columns):<2} | size={f.stat().st_size} bytes")

print("\n--- 3. Results Figures Inventory ---")
figs_dir = root / "results" / "figures"
for f in sorted(figs_dir.glob("*.png")):
    print(f"{f.name:<40} | size={f.stat().st_size} bytes")

print("\n--- 4. Source Code Files ---")
src_dir = root / "src" / "telotrade"
for f in sorted(src_dir.glob("*.py")):
    lines = len(f.read_text(encoding="utf-8").splitlines())
    print(f"{f.name:<25} | lines={lines:<5} | size={f.stat().st_size} bytes | hash={sha256_file(f)[:12]}...")
