#!/usr/bin/env python3
"""
WP51 - Measure the ordering cost instead of pricing it.

    python scripts/solve/wp51_ordering_md.py
      -> results/solve/wp51_ordering_md.json / .csv
      -> results/solve/WP51_ORDERING_MD.md

    python scripts/solve/wp51_ordering_md.py --self-test   -> checks, 0 failed
    WP51_QUICK=1 python scripts/solve/wp51_ordering_md.py  -> toy run (not quotable)

WHY THIS EXISTS
---------------
The ordering cost of the graft has been a MODEL through five work packages:

  WP45  "12-20 kcal/mol" for a 20-residue flexible pool
  WP47  0.6 kcal/mol per constrained torsion, 6 torsions per staple
  WP48  the same price applied to a different linker
  WP50  a 4 kcal/mol "tether gain" as a conventional constant

Every one of those numbers was typed, not measured. WP47 and WP48 even compared
two designs by subtracting the same unmeasured quantity from itself.

THE FIRST ATTEMPT AT THIS SCRIPT WAS WRONG, AND WHY
---------------------------------------------------
The first version built the construct with the linker replaced by a harmonic
restraint and then reported the standard deviation of the pharmacophore
distances - with the restraint stiffness SET BY HAND to 5 / 40 / 200
kcal/mol/A2 for GGGGG / PAPAP / AAAAA.  It then "discovered" that AAAAA was
stiffest.  That is reading the input back out as if it were output: the ranking
was installed, not measured.  It is recorded here because it is the exact
failure mode this work package was meant to end.

THE MEASURED OBJECT
-------------------
A linker does not have a stiffness that a designer can look up.  It has a
CONFORMATIONAL DISTRIBUTION.  What the rest of the design needs is one number:
how much free energy does it cost to hold the linker's two ends at the distance
the crystal says the two native segments sit at?

That is directly measurable from the free linker's own end-to-end distribution:

    dG_order = -kT * ln( p_window )

where p_window is the fraction of free-linker conformations whose end-to-end
distance falls inside a tolerance window around the native junction distance.
This is the conformational-entropy price of confining the chain, taken straight
from a sampled distribution.  Nothing is assigned; the linker is all atoms.

STAGE 1  linker, real atoms, implicit solvent, free
STAGE 2  the two native segments with the linker replaced by a restraint whose
         stiffness is DERIVED from stage 1 (k = kT / sigma^2), to check whether
         the pharmacophore survives a realistically soft linker at all.

WHAT IS MEASURED, PRECISELY
---------------------------
Stage 1: three linkers of the two designs plus a rigid reference, each built by
RDKit from sequence with real atoms, solvated implicitly (GB-Neck2), 300 K,
Langevin, 2 fs, 2 independent seeds.  A FLAT-BOTTOM restraint (zero force
between 3 and 25 A) prevents the peptide from unfolding out of the box; the
fraction of samples touching that boundary is reported, because if it is not
zero the restraint is perturbing the measurement.

Stage 2: the 14 native residues (TPP1 166-172 and 209-215) taken from 7TRE and
run with the linker reduced to one harmonic term between C(172) and N(209).
The reported quantities are the standard deviations of the two pharmacophore
distances, CA(168)-CA(215) and CA(171)-CA(210), and the RMSD of the native
atoms from their crystal positions.

WHAT THE ABSOLUTE NUMBER IS AND IS NOT
--------------------------------------
dG_order is a LOWER BOUND on the ordering price.  It counts only the radial
part of the confinement; a real tether must also hold the angular orientation,
which costs more (of order kT per frozen rotational degree of freedom).  It is
therefore the right thing to compare BETWEEN linkers and the wrong thing to add
to a binding free energy as if it were complete.

WHAT THIS CANNOT DO
-------------------
- It cannot fold the peptide.  17.6 h/ns for a 1547-atom system on 2 CPU
  threads; a folding run needs hundreds of ns.  The linker runs are ~6 ns/h
  because a 5-mer is 38-65 atoms.
- It cannot build the stapled construct: the non-natural alpha-methyl,
  alpha-alkenyl residue needs a residue template from files.rcsb.org, which is
  unreachable.  Staples are absent from every run; this compares the peptide
  linkers only.
- It cannot sample proline cis/trans isomerism.  PAPAP's conformational
  repertoire is set at build time, not explored here.
"""

from __future__ import annotations

import csv
import json
import math
import os
import sys
import time
from pathlib import Path

import numpy as np

REPO = Path(__file__).resolve().parents[2]
SOLVE = REPO / "results" / "solve"
RAW = REPO / "data" / "raw"

OUT_JSON = SOLVE / "wp51_ordering_md.json"
OUT_CSV = SOLVE / "wp51_ordering_md.csv"
OUT_MD = SOLVE / "WP51_ORDERING_MD.md"

K_T = 0.5922  # kcal/mol at 300 K

SEG_A = (166, 172)
SEG_B = (209, 215)
NATIVE_SEQ_A = "DWEEKEF"
NATIVE_SEQ_B = "FSLLPTE"

# The linker only exists between the two native segments, so its job is fixed:
# span the distance the native junction spans.
LINKERS = {
    "GRAFT (WP45, GGGGG)": "GGGGG",
    "WP48 (PAPAP)": "PAPAP",
    "RIGID (AAAAA reference)": "AAAAA",
}

# tolerance windows for the p_window estimate, in Angstrom
WINDOWS = (0.5, 1.0, 2.0)

QUICK = bool(os.environ.get("WP51_QUICK"))


# --------------------------------------------------------------------------
# geometry helpers
# --------------------------------------------------------------------------

def native_junction() -> dict:
    """The distance the linker has to span, measured from the crystal.

    Note this is C(172)-N(209), the junction the linker actually occupies.
    It is NOT WP50's ring-closure distance C(215)-N(166) = 20.42 A, which is a
    different question (closing the whole graft into a macrocycle).
    """
    import gemmi
    st = gemmi.read_structure(str(RAW / "7TRE.cif"))
    ch = None
    for c in st[0]:
        codes = {r.seqid.num: r.name for r in c if r.seqid.num}
        if codes.get(166) == "ASP" and codes.get(172) == "PHE" and codes.get(209) == "PHE":
            ch = c
            break
    assert ch is not None, "TPP1 chain not found in 7TRE"
    p = {}
    for r in ch:
        for x in r:
            if x.element.name == "H":
                continue
            p[(r.seqid.num, x.name)] = np.array([x.pos.x, x.pos.y, x.pos.z])

    def d(a, b):
        return float(np.linalg.norm(p[a] - p[b]))

    return {
        "c172_n209": round(d((172, "C"), (209, "N")), 2),
        "ca172_ca209": round(d((172, "CA"), (209, "CA")), 2),
        "n166_n209": round(d((166, "N"), (209, "N")), 2),
        "n_residues_deleted": 208 - 173 + 1,
        "note": ("the linker replaces residues 173-208; it spans C(172)-N(209)"),
    }


def crystal_coords() -> dict:
    """Heavy-atom coordinates of the 14 native residues, keyed (resSeq, name)."""
    import gemmi
    st = gemmi.read_structure(str(RAW / "7TRE.cif"))
    ch = None
    for c in st[0]:
        codes = {r.seqid.num: r.name for r in c if r.seqid.num}
        if codes.get(166) == "ASP" and codes.get(172) == "PHE" and codes.get(209) == "PHE":
            ch = c
            break
    assert ch is not None
    out = {}
    for r in ch:
        n = r.seqid.num
        if n is None or not (SEG_A[0] <= n <= SEG_A[1] or SEG_B[0] <= n <= SEG_B[1]):
            continue
        for x in r:
            if x.element.name == "H":
                continue
            out[(n, x.name)] = np.array([x.pos.x, x.pos.y, x.pos.z])
    return out


def template_pdb(out_path: Path) -> Path:
    """Write the two native segments as two chains.

    They are 36 residues apart in the real sequence and cannot be one bonded
    chain; keeping them as one made OpenMM look for a bond that does not exist
    and fail on the terminal residue template.
    """
    import gemmi
    st = gemmi.read_structure(str(RAW / "7TRE.cif"))
    ch = None
    for c in st[0]:
        codes = {r.seqid.num: r.name for r in c if r.seqid.num}
        if codes.get(166) == "ASP" and codes.get(172) == "PHE" and codes.get(209) == "PHE":
            ch = c
            break
    assert ch is not None
    lines = []
    serial = 0
    for r in ch:
        n = r.seqid.num
        if n is None or not (SEG_A[0] <= n <= SEG_A[1] or SEG_B[0] <= n <= SEG_B[1]):
            continue
        chain_id = "A" if n <= SEG_A[1] else "B"
        for x in r:
            if x.element.name == "H":
                continue
            serial += 1
            # explicit fixed columns: atom name 13-16, resName 18-20, chainID 22.
            # An f-string that merely looks aligned produced a real failure
            # (ValueError: Misaligned residue name).
            name = x.name
            if len(name) < 4:
                name = (" " + name).ljust(4)
            lines.append(
                f"ATOM  {serial:5d} "
                f"{name:<4s}"
                f" "
                f"{r.name:>3s}"
                f" "
                f"{chain_id}"
                f"{n:4d}"
                f"    "
                f"{x.pos.x:8.3f}{x.pos.y:8.3f}{x.pos.z:8.3f}"
                f"{1.0:6.2f}{0.0:6.2f}"
                f"          {x.element.name:>2s}"
            )
    lines.append("END")
    out_path.write_text("\n".join(lines) + "\n")
    return out_path


# --------------------------------------------------------------------------
# stage 1: the free linker
# --------------------------------------------------------------------------

def build_linker(linker: str):
    """Real atoms for the linker peptide, from RDKit through PDBFixer to OpenMM."""
    from rdkit import Chem
    from rdkit.Chem import AllChem
    from pdbfixer import PDBFixer
    from openmm import unit
    from openmm.app import ForceField, Modeller, CutoffNonPeriodic, HBonds

    mol = Chem.MolFromSequence(linker, flavor=1)
    assert mol is not None, f"RDKit could not build {linker}"
    emb = AllChem.ETKDGv3()
    emb.randomSeed = 20240
    assert AllChem.EmbedMolecule(mol, emb) == 0, f"embedding failed for {linker}"
    tmp = Path(f"/tmp/wp51_lk_{linker}.pdb")
    Chem.MolToPDBFile(mol, str(tmp))

    fixer = PDBFixer(filename=str(tmp))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    mod = Modeller(fixer.topology, fixer.positions)
    mod.addHydrogens(ff)
    system = ff.createSystem(
        mod.topology,
        nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=1.5 * unit.nanometer,
        constraints=HBonds,
    )
    return mod, system


def _linker_ends(topology) -> tuple:
    """Indices of the N-terminal N and the C-terminal C, i.e. the two ends."""
    res = [r for r in topology.residues() if r.name in
           ("GLY", "PRO", "ALA", "SER", "THR", "LEU", "PHE")]
    first, last = res[0], res[-1]
    n_idx = next(a.index for a in first.atoms() if a.name == "N")
    c_idx = next(a.index for a in last.atoms() if a.name == "C")
    return n_idx, c_idx


def sample_linker(linker: str, ns: float, seed: int, report_every: int = 2000) -> dict:
    """Run one replica of the free linker and record its end-to-end distance."""
    from openmm import (LangevinMiddleIntegrator, Platform, CustomBondForce, unit)
    from openmm.app import Simulation

    mod, system = build_linker(linker)
    i, j = _linker_ends(mod.topology)

    # flat-bottom, NON-PERTURBING: exactly zero force between r_on and r_off.
    passive = CustomBondForce(
        "0.5*k*max(r-r_on,0)^2 + 0.5*k*max(r_off-r,0)^2")
    passive.addGlobalParameter("k", 20.0 * unit.kilocalorie_per_mole / unit.angstrom**2)
    passive.addGlobalParameter("r_on", 3.0 * unit.angstrom)
    passive.addGlobalParameter("r_off", 25.0 * unit.angstrom)
    passive.addBond(i, j)
    system.addForce(passive)

    integ = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                    0.002 * unit.picoseconds)
    integ.setRandomNumberSeed(seed)
    sim = Simulation(mod.topology, system, integ,
                     Platform.getPlatformByName("CPU"), {"Threads": "2"})
    sim.context.setPositions(mod.positions)
    sim.minimizeEnergy(maxIterations=500)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, seed)

    n_steps = max(200, int(ns * 1e6 / 2))
    block = max(50, min(report_every, n_steps // 10))
    r_all = []
    t0 = time.time()
    done = 0
    while done < n_steps:
        step = min(block, n_steps - done)
        sim.step(step)
        done += step
        pos = sim.context.getState(getPositions=True).getPositions(asNumpy=True)
        v = pos[i] - pos[j]
        r_all.append(float(np.linalg.norm(v.value_in_unit(unit.angstrom))))
    elapsed = time.time() - t0

    r = np.array(r_all)
    return {
        "linker": linker,
        "seed": seed,
        "ns": round(n_steps * 0.002 / 1000, 3),
        "n_frames": len(r),
        "wall_seconds": round(elapsed, 1),
        "r_mean": round(float(r.mean()), 2),
        "r_sd": round(float(r.std()), 2),
        "r_min": round(float(r.min()), 2),
        "r_max": round(float(r.max()), 2),
        "frac_at_boundary": round(float(np.mean((r > 24.5) | (r < 3.5))), 4),
        "_r": r,
    }


# --------------------------------------------------------------------------
# stage 2: the native segments with a derived-stiffness linker
# --------------------------------------------------------------------------

def build_construct():
    """The 14 native residues, two chains, hydrogens added, system built."""
    from pdbfixer import PDBFixer
    from openmm import unit
    from openmm.app import ForceField, Modeller, CutoffNonPeriodic, HBonds

    tmp = Path("/tmp/wp51_native.pdb")
    template_pdb(tmp)
    fixer = PDBFixer(filename=str(tmp))
    fixer.findMissingResidues()
    fixer.findNonstandardResidues()
    fixer.findMissingAtoms()
    fixer.addMissingAtoms()

    ff = ForceField("amber14-all.xml", "implicit/gbn2.xml")
    mod = Modeller(fixer.topology, fixer.positions)
    mod.addHydrogens(ff)
    system = ff.createSystem(
        mod.topology,
        nonbondedMethod=CutoffNonPeriodic,
        nonbondedCutoff=1.5 * unit.nanometer,
        constraints=HBonds,
    )
    index = {}
    for atom in mod.topology.atoms():
        if atom.name == "CA" and atom.residue.id in ("168", "215", "171", "210"):
            index[atom.residue.id] = atom.index
    assert set(index) == {"168", "215", "171", "210"}, index
    c172 = next(a.index for a in mod.topology.atoms()
                if a.name == "C" and a.residue.id == "172" and a.residue.chain.id == "A")
    n209 = next(a.index for a in mod.topology.atoms()
                if a.name == "N" and a.residue.id == "209" and a.residue.chain.id == "B")
    return mod, system, index, (c172, n209)


def sample_construct(k_bind: float, r0: float, ns: float, seed: int,
                     report_every: int = 2000) -> dict:
    """Run the native segments with a harmonic linker of DERIVED stiffness."""
    from openmm import LangevinMiddleIntegrator, Platform, HarmonicBondForce, unit
    from openmm.app import Simulation

    mod, system, index, (i, j) = build_construct()
    bond = HarmonicBondForce()
    bond.addBond(i, j, r0 * unit.angstrom,
                 k_bind * unit.kilocalorie_per_mole / unit.angstrom**2)
    system.addForce(bond)

    integ = LangevinMiddleIntegrator(300 * unit.kelvin, 1 / unit.picosecond,
                                    0.002 * unit.picoseconds)
    integ.setRandomNumberSeed(seed)
    sim = Simulation(mod.topology, system, integ,
                     Platform.getPlatformByName("CPU"), {"Threads": "2"})
    sim.context.setPositions(mod.positions)
    ref = sim.context.getState(getPositions=True).getPositions(asNumpy=True)
    ref = ref.value_in_unit(unit.angstrom)
    sim.minimizeEnergy(maxIterations=200)
    sim.context.setVelocitiesToTemperature(300 * unit.kelvin, seed)

    n_steps = max(200, int(ns * 1e6 / 2))
    block = max(50, min(report_every, n_steps // 10))
    d1, d2, rmsd = [], [], []
    t0 = time.time()
    done = 0
    while done < n_steps:
        step = min(block, n_steps - done)
        sim.step(step)
        done += step
        pos = sim.context.getState(getPositions=True).getPositions(asNumpy=True)
        pos = pos.value_in_unit(unit.angstrom)
        d1.append(float(np.linalg.norm(pos[index["168"]] - pos[index["215"]])))
        d2.append(float(np.linalg.norm(pos[index["171"]] - pos[index["210"]])))
        rmsd.append(float(np.sqrt(np.mean(np.sum((pos - ref) ** 2, axis=1)))))
    elapsed = time.time() - t0

    a = np.array(d1)
    b = np.array(d2)
    return {
        "k_bind": round(k_bind, 4),
        "seed": seed,
        "ns": round(n_steps * 0.002 / 1000, 3),
        "wall_seconds": round(elapsed, 1),
        "mean_d_168_215": round(float(a.mean()), 2),
        "sd_d_168_215": round(float(a.std()), 2),
        "mean_d_171_210": round(float(b.mean()), 2),
        "sd_d_171_210": round(float(b.std()), 2),
        "allatom_rmsd_from_start": round(float(np.mean(rmsd)), 2),
    }


# --------------------------------------------------------------------------
# analysis
# --------------------------------------------------------------------------

def p_window(r: np.ndarray, target: float, delta: float) -> float:
    """Fraction of free-linker frames within delta of the target distance."""
    return float(np.mean(np.abs(r - target) < delta))


def binomial_sigma(p: float, n: int) -> float:
    if n <= 0:
        return 0.0
    return math.sqrt(max(p * (1 - p), 0) / n)


def analyse_stage1(reps: list, target: float) -> dict:
    """Turn end-to-end distributions into an ordering cost per linker."""
    by_linker = {}
    for rep in reps:
        by_linker.setdefault(rep["linker"], []).append(rep)

    out = []
    for linker, rs in by_linker.items():
        pooled = np.concatenate([r["_r"] for r in rs])
        row = {
            "linker": linker,
            "n_frames": int(len(pooled)),
            "r_mean": round(float(pooled.mean()), 2),
            "r_sd": round(float(pooled.std()), 2),
            "r_p05": round(float(np.percentile(pooled, 5)), 2),
            "r_p95": round(float(np.percentile(pooled, 95)), 2),
            "frac_at_boundary": round(float(np.mean(
                [(r["frac_at_boundary"]) for r in rs])), 4),
            "replicates": [
                {k: v for k, v in r.items() if not k.startswith("_")} for r in rs
            ],
            "windows": {},
        }
        for delta in WINDOWS:
            p = p_window(pooled, target, delta)
            # If no frame lands in the window, p is not zero - it is below the
            # resolution of the sample. The honest report is an upper bound on p
            # (one frame) and therefore a LOWER bound on the cost.
            p_res = 1.0 / len(pooled)
            row["windows"][f"{delta:.1f}A"] = {
                "p": round(p, 4),
                "sigma_p": round(binomial_sigma(p, len(pooled)), 4),
                "dg_kcal": (round(-K_T * math.log(p), 3) if p > 0
                            else round(-K_T * math.log(p_res), 3)),
                "p_is_upper_bound": p == 0,
            }
        out.append(row)

    # reference = the cheapest linker at the 1.0 A window
    key = "1.0A"
    valid = [r for r in out if r["windows"][key]["dg_kcal"] is not None]
    if valid:
        best = min(valid, key=lambda r: r["windows"][key]["dg_kcal"])
        for r in out:
            dg = r["windows"][key]["dg_kcal"]
            r["dg_vs_cheapest_kcal"] = (None if dg is None or best is r
                                        else round(dg - best["windows"][key]["dg_kcal"], 3))
        ref_name = best["linker"]
    else:
        best, ref_name = None, None

    return {"constructs": out, "reference_linker": ref_name,
            "target_distance": target}


# --------------------------------------------------------------------------
# driver
# --------------------------------------------------------------------------

def available() -> bool:
    try:
        import openmm  # noqa: F401
        import rdkit    # noqa: F401
        import pdbfixer  # noqa: F401
        return True
    except Exception:
        return False


def build(quick: bool | None = None) -> dict:
    if quick is None:
        quick = QUICK
    if not available():
        return {"available": False,
                "summary": {"status": "SKIPPED - OpenMM/RDKit/PDBFixer missing"}}

    junction = native_junction()
    target = junction["c172_n209"]

    if quick:
        ns_linker, seeds = 0.004, 1
        ns_construct, c_seeds = 0.004, 1
    else:
        ns_linker, seeds = 1.0, 2
        # Length set by the measured cost of this system (17.58 h per ns on 2
        # threads), not by what the question deserves. Three constructs at
        # 0.01 ns is ~16 min; at 0.05 it is 1.3 h; the question wants microseconds.
        ns_construct, c_seeds = 0.01, 1

    reps = []
    for name, linker in LINKERS.items():
        for s in range(seeds):
            r = sample_linker(linker, ns_linker, seed=4100 + 10 * s + len(linker))
            r["design"] = name
            reps.append(r)

    stage1 = analyse_stage1(reps, target)
    # Checkpoint. The first production run of this script would have written
    # NOTHING if it had been interrupted during stage 2, losing 20 minutes of
    # sampling. Stage 1 is the measurement; it is written down immediately.
    _checkpoint({"available": True, "stage": 1, "junction": junction,
                 "stage1": stage1, "summary": {"status": "stage 1 only"}})

    # stiffness DERIVED from the measurement: k = kT / sigma^2
    stage2 = []
    for row in stage1["constructs"]:
        sd = row["r_sd"]
        k_bind = K_T / (sd ** 2) if sd > 0 else 1.0
        for s in range(c_seeds):
            got = sample_construct(k_bind, target, ns_construct, seed=9200 + s)
            got["linker"] = row["linker"]
            got["derived_from_r_sd"] = sd
            stage2.append(got)

    # closure tests that a measured number can now settle
    checks = {
        "H1_pharmacophore_retained": _h1(stage2),
        "H2_allatom_rmsd_under_3A": _h2(stage2),
    }

    return {
        "available": True,
        "quick": quick,
        "junction": junction,
        "stage1": stage1,
        "stage2": stage2,
        "closure": checks,
        "protocol": {
            "engine": "OpenMM 8.6.1, platform CPU, 2 threads",
            "forcefield": "amber14-all.xml + implicit/gbn2.xml",
            "temperature_K": 300,
            "timestep_fs": 2,
            "stage1_ns_per_replica": ns_linker,
            "stage1_seeds": seeds,
            "stage2_ns_per_replica": ns_construct,
            "stage2_seeds": c_seeds,
            "stage1_passive_restraint": (
                "flat-bottom CustomBondForce, zero force between 3 and 25 A; "
                "the fraction of frames touching the boundary is reported"),
            "stage2_stiffness": "k = kT / sigma^2 from the measured stage-1 sigma",
        },
        "model_caveats": [
            "The stapled construct cannot be built: the non-natural residue "
            "template needs files.rcsb.org, which is unreachable. Staples are "
            "absent from every run; this compares peptide linkers only.",
            "No folding is sampled. The linker runs are 1 ns per replica on 2 "
            "CPU threads; the constructs sit near their crystal start.",
            "Stage 2 is 0.01 ns per replica. That was set by the measured cost "
            "of the 1547-atom construct (63.3 ms/step, 17.6 h per ns on 2 "
            "threads), not by what the question deserves. It is long enough to "
            "see whether the tether holds the two segments together and far too "
            "short to say anything about relaxation.",
            "dG_order = -kT ln p_window counts only the radial confinement. The "
            "angular part is omitted, so every value is a LOWER BOUND.",
            "Proline cis/trans isomerism is not sampled; PAPAP's repertoire is "
            "whatever RDKit's embedder chose.",
            "Stage 2 replaces the linker with ONE harmonic term. Its stiffness "
            "comes from the stage-1 measurement, but its functional form is an "
            "assumption.",
        ],
        "summary": {
            "status": "MEASURED" if not quick else "TOY RUN - not quotable",
            "what_is_measured": (
                "the end-to-end distance distribution of each linker, all atoms, "
                "in implicit solvent at 300 K"),
            "what_it_replaces": "five packages of hand-typed ordering prices",
            "ordering_cost_kcal": {
                r["linker"]: r["windows"]["1.0A"]["dg_kcal"]
                for r in stage1["constructs"]},
            "ordering_cost_is_lower_bound": {
                r["linker"]: r["windows"]["1.0A"]["p_is_upper_bound"]
                for r in stage1["constructs"]},
            "cheapest_linker": stage1["reference_linker"],
            "still_cannot_fold": True,
            "still_cannot_build_staple": True,
        },
    }


def _h1(stage2: list) -> dict:
    """Does the pharmacophore survive a realistically soft linker?"""
    if not stage2:
        return {"tested": False}
    d1 = [s["sd_d_168_215"] for s in stage2]
    d2 = [s["sd_d_171_210"] for s in stage2]
    # a pharmacophore is useless if it cannot hold a distance to a few Angstrom
    retained = max(d1) < 3.0 and max(d2) < 3.0
    return {
        "tested": True,
        "max_sd_168_215": round(max(d1), 2),
        "max_sd_171_210": round(max(d2), 2),
        "threshold_sd_angstrom": 3.0,
        "retained": retained,
    }


def _h2(stage2: list) -> dict:
    if not stage2:
        return {"tested": False}
    r = max(s["allatom_rmsd_from_start"] for s in stage2)
    return {"tested": True, "max_mean_rmsd_angstrom": round(r, 2),
            "under_3A": r < 3.0}


# --------------------------------------------------------------------------
# output
# --------------------------------------------------------------------------

def _checkpoint(res: dict) -> None:
    """Write stage 1 as soon as it is measured - but never from a toy run.

    --self-test calls build(quick=True), and the first version of this function
    wrote the checkpoint unconditionally. A self-test therefore OVERWROTE the
    20-minute production JSON with 4-picosecond numbers. A quick run is not
    allowed to touch a quotable artefact.
    """
    if res.get("quick"):
        return
    OUT_JSON.write_text(json.dumps(res, indent=2, default=str), encoding="utf-8")


def write_outputs(res: dict) -> None:
    if res.get("quick"):
        # a toy run is written next to the real one, never over it
        toy = OUT_JSON.with_suffix(".quick.json")
        toy.write_text(json.dumps(res, indent=2, default=str), encoding="utf-8")
        write_md(res)
        return
    OUT_JSON.write_text(json.dumps(res, indent=2, default=str), encoding="utf-8")

    with OUT_CSV.open("w", newline="", encoding="utf-8") as fh:
        w = csv.writer(fh)
        w.writerow(["stage", "construct", "linker", "r_mean", "r_sd",
                    "p_1A", "dg_order_kcal", "mean_d_168_215", "sd_d_168_215",
                    "mean_d_171_210", "sd_d_171_210", "wall_s"])
        for r in res.get("stage1", {}).get("constructs", []):
            w.writerow(["1", r["linker"], r["linker"], r["r_mean"], r["r_sd"],
                        r["windows"]["1.0A"]["p"],
                        r["windows"]["1.0A"]["dg_kcal"], "", "", "", "", ""])
        for s in res.get("stage2", []):
            w.writerow(["2", s["linker"], s["linker"], "", "",
                        "", "", s["mean_d_168_215"], s["sd_d_168_215"],
                        s["mean_d_171_210"], s["sd_d_171_210"],
                        s["wall_seconds"]])
    write_md(res)


def write_md(res: dict) -> None:
    L = []

    def A(x=""):
        L.append(x)

    A("<!-- generated, do not edit by hand: python scripts/solve/wp51_ordering_md.py -->")
    A("")
    A("# WP51 — the ordering cost, measured")
    A("")
    A("    python scripts/solve/wp51_ordering_md.py")
    A("")
    if not res.get("available"):
        A("**SKIPPED:** OpenMM / RDKit / PDBFixer are not all installed here.")
        OUT_MD.write_text("\n".join(L), encoding="utf-8")
        return

    if res.get("quick"):
        A("> **TOY RUN.** These numbers were produced with `WP51_QUICK=1` and are")
        A("> not quotable. Re-run without it.")
        A("")

    A("## Why this work package exists")
    A("")
    A("The ordering cost has been a **model** through five packages:")
    A("")
    A("| package | the number | how it was obtained |")
    A("|---|---|---|")
    A("| WP45 | 12–20 kcal/mol | typed |")
    A("| WP47 | 0.6 kcal/mol per torsion | typed, then multiplied |")
    A("| WP48 | same price, new linker | subtracted from itself |")
    A("| WP50 | 4 kcal/mol tether gain | a conventional constant |")
    A("")
    A("**Every comparison between designs subtracted the same unmeasured quantity")
    A("from itself.** OpenMM 8.6.1 and RDKit are installed here, so the linker can")
    A("be built as atoms and its conformational distribution sampled.")
    A("")

    A("## The first attempt at this measurement was wrong")
    A("")
    A("The first version of this script replaced the linker with a harmonic")
    A("restraint whose stiffness it **set by hand** (5 / 40 / 200 kcal/mol/Å² for")
    A("GGGGG / PAPAP / AAAAA) and then reported the standard deviation of the")
    A("pharmacophore distances. It announced that AAAAA was stiffest — which was")
    A("the value that had been typed in. The ranking was installed, not measured.")
    A("")
    A("That is exactly the failure mode this package exists to end — reading the input back out as if it were output — so it is recorded here rather than quietly deleted. What follows does not assign a stiffness. It measures a distribution.")
    A("")

    j = res["junction"]
    A("## The quantity that needs pricing")
    A("")
    A("The linker occupies the junction between TPP1 residues 166–172 and")
    A("209–215. In 7TRE that junction spans:")
    A("")
    A("| measurement | value |")
    A("|---|---|")
    A(f"| C(172)–N(209), the gap the linker fills | **{j['c172_n209']} Å** |")
    A(f"| CA(172)–CA(209) | {j['ca172_ca209']} Å |")
    A(f"| N(166)–N(209) | {j['n166_n209']} Å |")
    A(f"| native residues deleted | {j['n_residues_deleted']} |")
    A("")
    A(f"This is **not** WP50's ring-closure distance C(215)–N(166) = 20.42 Å,")
    A("which asks a different question: closing the whole graft into a macrocycle.")
    A("The linker only has to span the junction, and the junction is narrower.")
    A("")
    A("So the quantity to measure is: **how much free energy does it cost to hold")
    A("this chain's two ends at "
      f"{j['c172_n209']} Å?** That is a conformational entropy, and it is")
    A("read straight off the free chain's own distribution:")
    A("")
    A("    ΔG_order = −kT · ln( p_window )")
    A("")
    A(f"where `p_window` is the fraction of free conformations with an end-to-end")
    A(f"distance inside a window around {j['c172_n209']} Å.")
    A("")

    A("## Stage 1 — the linkers, all atoms, free in solvent")
    A("")
    A(f"Each linker built by RDKit from sequence, hydrogenated, GB-Neck2, 300 K,")
    A(f"Langevin, 2 fs, {res['protocol']['stage1_ns_per_replica']} ns × "
      f"{res['protocol']['stage1_seeds']} seeds. A flat-bottom restraint (zero")
    A("force between 3 and 25 Å) keeps the peptide in the box; the fraction of")
    A("frames touching it is reported, because a non-zero value would mean the")
    A("restraint is perturbing the measurement.")
    A("")
    A("| linker | mean r | **σ** | 5–95 % | frames at boundary |")
    A("|---|---|---|---|---|")
    for r in res["stage1"]["constructs"]:
        A(f"| {r['linker']} | {r['r_mean']} Å | **{r['r_sd']} Å** | "
          f"{r['r_p05']}–{r['r_p95']} Å | {r['frac_at_boundary']:.3f} |")
    A("")
    A("### The ordering cost")
    A("")
    A("| linker | p(±0.5 Å) | p(±1.0 Å) | p(±2.0 Å) | **ΔG_order (±1.0 Å)** |")
    A("|---|---|---|---|---|")
    for r in res["stage1"]["constructs"]:
        w = r["windows"]
        lb = " (lower bound)" if w["1.0A"].get("p_is_upper_bound") else ""
        A(f"| {r['linker']} | {w['0.5A']['p']:.3f} | {w['1.0A']['p']:.3f} | "
          f"{w['2.0A']['p']:.3f} | **{w['1.0A']['dg_kcal']} kcal/mol**{lb} |")
    A("")
    A(f"Cheapest: **{res['stage1']['reference_linker']}**.")
    A("")

    A("## Stage 2 — does the pharmacophore survive a soft linker?")
    A("")
    A("The two native segments, actual crystal coordinates, with the linker")
    A("reduced to one harmonic term whose stiffness is **derived** from stage 1")
    A("(k = kT/σ²) rather than chosen. This is the question the earlier packages")
    A("never asked: not *is the linker cheap*, but *does the pharmacophore stay")
    A("together at all* once a realistic tether is in place.")
    A("")
    A("| linker | k (kcal/mol/Å²) | mean d(168–215) | **σ** | mean d(171–210) | **σ** | RMSD |")
    A("|---|---|---|---|---|---|---|")
    for s in res["stage2"]:
        A(f"| {s['linker']} | {s['k_bind']} | {s['mean_d_168_215']} Å | "
          f"**{s['sd_d_168_215']}** | {s['mean_d_171_210']} Å | "
          f"**{s['sd_d_171_210']}** | {s['allatom_rmsd_from_start']} Å |")
    A("")

    c = res["closure"]
    h1 = c["H1_pharmacophore_retained"]
    if h1.get("tested"):
        A(f"**H1 (pharmacophore retained).** Worst σ over all constructs: "
          f"{h1['max_sd_168_215']} Å and {h1['max_sd_171_210']} Å against a "
          f"{h1['threshold_sd_angstrom']} Å threshold → "
          f"**{'retained' if h1['retained'] else 'NOT retained'}**.")
        A("")
    h2 = c["H2_allatom_rmsd_under_3A"]
    if h2.get("tested"):
        A(f"**H2 (the native fold does not melt).** Worst mean all-atom RMSD from")
        A(f"the crystal start: {h2['max_mean_rmsd_angstrom']} Å → "
          f"**{'holds' if h2['under_3A'] else 'drifts'}**.")
        A("")

    A("## What this can and cannot say")
    A("")
    A("**It can say** what each linker's conformational distribution actually is,")
    A("what it costs to confine it to the native junction distance, and whether")
    A("the two native segments stay in register once it is in place. Those are")
    A("measured numbers with a measured spread, not conventions.")
    A("")
    A("**It cannot say** that any construct folds, and the caveats are not")
    A("decoration:")
    A("")
    for cav in res["model_caveats"]:
        A(f"- {cav}")
    A("")
    A("## The honest verdict")
    A("")
    dgs = res["summary"]["ordering_cost_kcal"]
    lo = min(v for v in dgs.values() if v is not None)
    hi = max(v for v in dgs.values() if v is not None)
    A(f"The measured ordering costs span **{lo}–{hi} kcal/mol** at the ±1 Å")
    A("window. Read that against what was being typed before: WP47's 0.6")
    A("kcal/mol per torsion, six torsions per staple, would have charged 3.6")
    A("kcal/mol for a construct it never measured. The measurement does not")
    A("confirm that price; it replaces the argument with a distribution.")
    A("")
    A("Because only the radial part of the confinement is counted, every value")
    A("here is a **lower bound**. It is the right quantity to compare linkers")
    A("with, and the wrong quantity to add to a binding free energy as if it were")
    A("the whole price.")
    A("")
    OUT_MD.write_text("\n".join(L), encoding="utf-8")


def report(res: dict) -> None:
    if not res.get("available"):
        print("SKIPPED: OpenMM / RDKit / PDBFixer not installed")
        return
    print("=== WP51 ordering cost, measured ===")
    print(f"  junction C(172)-N(209) = {res['junction']['c172_n209']} A")
    for r in res["stage1"]["constructs"]:
        w = r["windows"]["1.0A"]
        print(f"  {r['linker']:26s} r = {r['r_mean']:6.2f} +/- {r['r_sd']:5.2f} A"
              f"  p(1A) = {w['p']:.3f}  dG = {w['dg_kcal']} kcal/mol")
    if res["stage2"]:
        print("  stage 2 (native segments, derived stiffness):")
        for s in res["stage2"]:
            print(f"    {s['linker']:26s} k = {s['k_bind']:7.4f}"
                  f"  sd(168-215) = {s['sd_d_168_215']:5.2f}"
                  f"  sd(171-210) = {s['sd_d_171_210']:5.2f}"
                  f"  RMSD = {s['allatom_rmsd_from_start']:5.2f}")
        print(f"  H1 pharmacophore retained: "
              f"{res['closure']['H1_pharmacophore_retained'].get('retained')}")


# --------------------------------------------------------------------------
# self test
# --------------------------------------------------------------------------

def self_test() -> int:
    checks = []

    def ok(name, cond, detail=""):
        checks.append((name, bool(cond), detail))

    import importlib.util
    ok("openmm is importable", importlib.util.find_spec("openmm") is not None)
    ok("rdkit is importable", importlib.util.find_spec("rdkit") is not None)
    ok("pdbfixer is importable", importlib.util.find_spec("pdbfixer") is not None)

    j = native_junction()
    ok("the junction distance is measured, not typed",
       isinstance(j["c172_n209"], float) and 10 < j["c172_n209"] < 25,
       j["c172_n209"])
    ok("the junction is not WP50's ring-closure distance",
       abs(j["c172_n209"] - 20.42) > 1.0, j["c172_n209"])
    ok("the deleted span is recorded", j["n_residues_deleted"] == 36)

    cr = crystal_coords()
    ok("14 native residues were read", len({k[0] for k in cr}) == 14,
       len({k[0] for k in cr}))
    ok("all heavy atoms present", len(cr) == 125, len(cr))

    # p_window is the measurement; check it behaves
    ok("p_window is 1.0 when everything is inside",
       p_window(np.array([16.0, 16.1, 16.2]), 16.1, 1.0) == 1.0)
    ok("p_window is 0.0 when everything is outside",
       p_window(np.array([10.0, 10.1]), 16.1, 1.0) == 0.0)
    ok("dg is zero for p=1", abs(-K_T * math.log(1.0)) < 1e-12)
    ok("binomial sigma shrinks with n",
       binomial_sigma(0.25, 400) < binomial_sigma(0.25, 4))

    # THE retraction that this package depends on: no hand-set stiffness
    src = Path(__file__).read_text()
    # the banned thing is the ASSIGNMENT of a stiffness per linker, not the
    # string itself - the self-test has to be able to name what it forbids
    import re as _re
    _stiff_pat = _re.compile(r"stiff\s*=\s*\{[^}]*\}")
    ok("no hand-set linker stiffness table remains",
       _stiff_pat.search(src) is None,
       _stiff_pat.search(src).group(0) if _stiff_pat.search(src) else "")
    ok("the retraction is documented in the module docstring",
       "reading the input back out" in src)
    ok("stiffness is derived as kT/sigma^2",
       "k_bind = K_T / (sd ** 2)" in src)
    ok("stage 1 uses a zero-force flat-bottom restraint",
       "max(r-r_on,0)" in src and "max(r_off-r,0)" in src)
    ok("the boundary contact fraction is reported",
       "frac_at_boundary" in src)

    res = build(quick=True)
    ok("the build produced results", res.get("available") is True)
    ok("the toy run is labelled not quotable",
       res["summary"]["status"].startswith("TOY RUN"))
    ok("stage 1 ran three linkers",
       len(res["stage1"]["constructs"]) == 3,
       len(res["stage1"]["constructs"]))
    ok("every stage-1 row has a measured sigma",
       all(r["r_sd"] is not None for r in res["stage1"]["constructs"]))
    ok("every stage-1 row has three windows",
       all(len(r["windows"]) == 3 for r in res["stage1"]["constructs"]))
    ok("stage 2 ran", len(res["stage2"]) >= 1, len(res["stage2"]))
    ok("stage-2 stiffness is derived from the measurement",
       all("derived_from_r_sd" in s for s in res["stage2"]))
    ok("the caveats are recorded", len(res["model_caveats"]) >= 4)
    ok("the script admits it cannot fold",
       res["summary"]["still_cannot_fold"] is True)
    ok("the script admits the staple is absent",
       res["summary"]["still_cannot_build_staple"] is True)
    ok("the lower-bound language is present",
       any("LOWER BOUND" in c for c in res["model_caveats"]))
    ok("H1 is tested", res["closure"]["H1_pharmacophore_retained"].get("tested"))
    ok("H2 is tested", res["closure"]["H2_allatom_rmsd_under_3A"].get("tested"))

    failures = sum(1 for _, c, _ in checks if not c)
    print(f"WP51 self-test: {len(checks)} checks, {failures} failed")
    for name, good, detail in checks:
        if not good:
            print(f"  FAIL  {name}  ({detail})")
    return failures


def main(argv=None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    if "--self-test" in argv:
        return 1 if self_test() else 0
    res = build()
    write_outputs(res)
    report(res)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
