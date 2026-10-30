# Milestone 5 Remediation Explorer 3 Report: Resource Cleanup, Atomicity & Figures 1/3 Audit

**Explorer**: `explorer_m5_fix_3`  
**Parent Orchestrator**: `eb46fecb-7562-4ce7-aa98-998216cee2c4`  
**Date**: 2026-10-16  
**Target Codebase**: `src/telotrade/figures.py`  
**Objective**:
1. Audit Figure 1 and Figure 3 for covert hardcoding or synthetic fallback data.
2. Analyze resource leaks (`plt.close(fig)` vs parameterless `plt.close()`) and formulate exception-safe `try...finally: plt.close(fig)` design across all 4 figure generators.
3. Eliminate buffering issues and intermittent 0-byte file generation with atomic safe saving and buffer flushing (`os.fsync` + `os.replace`).
4. Provide concrete code patterns for Worker remediation.

---

## 1. Observation

### Observation 1: Covert Hardcoded Values in Figure 1 & Figure 3 Annotations
In `src/telotrade/figures.py`:
1. **Figure 1 (`generate_figure1_interface_contacts`) lines 131–163**:
   The annotation text strings hardcode specific measured distance numbers:
   ```python
   # Line 133
   ax1.annotate(
       "Arg774 -- Glu168\n(2.32 Å Salt Bridge)",
       xy=(168, 774),
       ...
   )
   # Line 144
   ax1.annotate(
       "Lys78 -- Glu215\n(3.58 Å Latch)",
       xy=(215, 78),
       ...
   )
   ```
   While `pair_summary` dynamically aggregates contacts from `df_contacts`, these specific text labels hardcode `2.32 Å` and `3.58 Å` rather than formatting the value dynamically from `pair_summary`.
2. **Figure 3 (`generate_figure3_tpp1_conformational_shifts`) lines 386–433**:
   Both target point coordinates (`xy=(res, dist)`) and annotation text strings are hardcoded with literal floating-point constants:
   ```python
   # Line 388-389
   ax1.annotate(
       "Lys170\n(Peak: 4.29 Å)",
       xy=(170, 4.293),
       ...
   )
   # Line 398-399
   ax1.annotate(
       "Glu169\n(3.81 Å)",
       xy=(169, 3.809),
       ...
   )
   # Line 407-408
   ax1.annotate(
       "Glu168\n(1.25 Å)",
       xy=(168, 1.253),
       ...
   )
   # Line 416-417
   ax1.annotate(
       "Leu183 Rigid Anchor\n(0.51 Å Displacement)",
       xy=(183, 0.512),
       ...
   )
   # Line 426-427
   ax1.annotate(
       "Leu212\n(1.04 Å)",
       xy=(212, 1.040),
       ...
   )
   ```
   If `tpp1_bound_vs_unbound_rmsd.csv` values are updated or perturbed, the line and scatter points shift, but the arrows point to the old hardcoded coordinates `(170, 4.293)` in empty space, and the callout text displays outdated numbers.
3. In line 450, `ax2.set_xlim(88, 243)` hardcodes fixed x-axis limits instead of deriving them from `resnums.min() - 2` and `resnums.max() + 2`.

### Observation 2: Exception-Unsafe Figure Lifecycle Across All 4 Generators
In `src/telotrade/figures.py`:
- Line 214: `plt.close()` in `generate_figure1_interface_contacts`
- Line 320: `plt.close()` in `generate_figure2_cross_state_rmsd`
- Line 454: `plt.close()` in `generate_figure3_tpp1_conformational_shifts`
- Line 578: `plt.close()` in `generate_figure4_cutoff_robustness`

Every figure generation function instantiates a figure via `plt.subplots()` or `plt.figure()`, but:
1. There is NO `try...finally` block. If any exception occurs during data filtering, plotting, or file saving, `plt.close()` is never executed. Matplotlib retains the figure in `matplotlib._pylab_helpers.Gcf.figs`, leaking memory across batch runs and test suites.
2. The call `plt.close()` without arguments targets `Gcf.get_active()`. If the active figure was modified by another thread or operation, the target figure handle `fig` remains unclosed.

### Observation 3: Intermittent 0-Byte Artifact Flakiness
In `tests/e2e/test_tier4_scenarios.py` line 134:
```python
for fig in figures:
    assert fig.stat().st_size > 1024, f"Figure {fig.name} is unexpectedly small"
```
During test runs, `figure4_cutoff_robustness.png` failed with `st_size == 0`.
Cause:
`plt.savefig(output_path)` opens the destination file in `"wb"` mode before rendering the canvas.
If an exception or early abort occurs after file creation, or if Windows filesystem write caching delays committing data to disk before `test_tier4_scenarios.py` inspects `st_size`, the file is observed as 0 bytes. Furthermore, `generate_all_figures` (lines 582–612) performs no verification of file existence or file size before returning.

---

## 2. Logic Chain

1. **Covert Hardcoding & Maintenance Risk**:
   - In Figure 1, the annotation text strings embed numbers `2.32 Å` and `3.58 Å`.
   - In Figure 3, both arrow target coordinates `xy=(res, dist)` and labels embed `4.293`, `3.809`, `1.253`, `0.512`, and `1.040`.
   - If upstream coordinates or alignments change, the plot points move while the annotations remain stuck on the hardcoded positions, breaking graphical alignment and reporting inaccurate numbers.
   - Replacing these with dynamic lookups (`pair_summary` in Figure 1, `df_tpp1.set_index("auth_seq_id")` in Figure 3) eliminates this desynchronization.

2. **Resource Management & Leak Prevention**:
   - Creating figures via `plt.subplots()` allocates memory in C++ / AGG backend buffers and registers the figure globally in `Gcf`.
   - Enclosing the entire plotting and file saving workflow in `try: ... finally: plt.close(fig)` guarantees that `fig` is cleanly destroyed, even if an exception is raised, eliminating memory leaks and `RuntimeWarning: More than 20 figures have been opened`.

3. **Atomic File Commitment & 0-Byte Elimination**:
   - Saving directly to `output_path` truncates the file immediately. If writing is incomplete or buffers are uncommitted, `stat().st_size` can read 0 bytes.
   - Introducing an atomic saving helper `save_figure_safely`:
     * Writes to a PID-tagged temporary file in the destination folder (`.{name}.{pid}.tmp`).
     * Flushes the Python file buffer (`fh.flush()`) and executes `os.fsync(fh.fileno())` to force the OS kernel to flush write caches to disk.
     * Validates that `tmp_path.stat().st_size > 1024` before committing.
     * Atomically moves the temporary file to `output_path` via `os.replace`. On Windows, `os.replace` replaces the destination atomically without ever leaving a 0-byte file in place.
     * Cleans up any dangling temporary files in a `finally` block.
   - Adding a post-generation integrity check in `generate_all_figures` asserts that all 4 figures exist and exceed 1024 bytes.

---

## 3. Caveats

1. In Figure 1, the vertical guide lines at `tel_patch_x = [168, 169, 170, 171, 183, 212]` represent the canonical TEL-patch definition from literature (`PROJECT.md` Feature 10). They can remain defined or be dynamically populated via `sorted(df_contacts[df_contacts["is_tel_patch"]]["tpp1_resnum"].unique())`.
2. In Figure 3, the classification bands (`ax1.axhspan(0, 1.0, ...)` etc.) and DE-loop highlight (`ax1.axvspan(167.5, 171.5, ...)`) are structural domain boundaries established by the scientific protocol, not data-dependent metrics.
3. This report focuses on resource cleanup, atomicity, and Figure 1/3 audits. Figures 2 and 4 data query implementations are also being addressed by companion explorers `explorer_m5_fix_1` and `explorer_m5_fix_2`.

---

## 4. Conclusion & Recommended Worker Implementation Plan

Worker should apply the following changes to `src/telotrade/figures.py`:

### Step 1: Add `import os` and Atomic Helper `save_figure_safely`
Place `save_figure_safely` right after `set_publication_style()`:
```python
import os
import time
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd


def save_figure_safely(fig: plt.Figure, output_path: Path | str, dpi: int = 300) -> Path:
    """
    Atomically save a matplotlib figure and guarantee non-zero file size (> 1024 bytes).

    Guarantees:
    1. Writes to a temporary file on the same filesystem/directory.
    2. Flushes Python stream buffers and forces OS disk cache commit via os.fsync.
    3. Validates file size > 1024 bytes before commitment.
    4. Atomically replaces target path via os.replace (preventing 0-byte visibility).
    5. Cleans up temporary file if any error occurs.
    """
    output_path = Path(output_path).resolve()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    tmp_path = output_path.with_name(f".{output_path.name}.{os.getpid()}_{time.time_ns()}.tmp")

    try:
        with open(tmp_path, "wb") as fh:
            fig.savefig(fh, format="png", dpi=dpi, bbox_inches="tight")
            fh.flush()
            os.fsync(fh.fileno())

        size = tmp_path.stat().st_size
        if size <= 1024:
            raise IOError(
                f"Generated figure {output_path.name} is unexpectedly small ({size} bytes, expected > 1024 bytes)"
            )

        os.replace(tmp_path, output_path)

        committed_size = output_path.stat().st_size
        if committed_size <= 1024:
            raise IOError(
                f"Committed figure {output_path.name} has invalid size ({committed_size} bytes)"
            )
    finally:
        if tmp_path.exists():
            try:
                tmp_path.unlink()
            except OSError:
                pass

    return output_path
```

### Step 2: Remediate Figure 1 (`generate_figure1_interface_contacts`)
1. Wrap the entire plotting body in `try: ... finally: plt.close(fig)`.
2. Dynamically extract annotation distances:
```python
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6.5), gridspec_kw={"width_ratios": [1.2, 1]})
    try:
        # Panel A plotting ...
        
        # Helper to query pair distance dynamically
        def _get_pair_dist(t_res: int, p_res: int) -> float | None:
            m = pair_summary[(pair_summary["tert_resnum"] == t_res) & (pair_summary["tpp1_resnum"] == p_res)]
            return float(m["min_dist"].min()) if not m.empty else None

        # Highlight Key Salt Bridge
        sb_dist = _get_pair_dist(774, 168)
        sb_text = f"Arg774 -- Glu168\n({sb_dist:.2f} Å Salt Bridge)" if sb_dist is not None else "Arg774 -- Glu168\n(Salt Bridge)"
        ax1.annotate(
            sb_text,
            xy=(168, 774),
            xytext=(145, 830),
            arrowprops={"facecolor": "#D90429", "shrink": 0.08, "width": 1.5, "headwidth": 6},
            fontsize=8.5,
            fontweight="bold",
            bbox={"boxstyle": "round,pad=0.3", "fc": "#FFF0F3", "ec": "#D90429", "lw": 1},
        )

        # Highlight Lys78 -- Glu215
        latch_dist = _get_pair_dist(78, 215)
        latch_text = f"Lys78 -- Glu215\n({latch_dist:.2f} Å Latch)" if latch_dist is not None else "Lys78 -- Glu215\n(Latch)"
        ax1.annotate(
            latch_text,
            xy=(215, 78),
            xytext=(185, 130),
            arrowprops={"facecolor": "#D90429", "shrink": 0.08, "width": 1.5, "headwidth": 6},
            fontsize=8.5,
            fontweight="bold",
            bbox={"boxstyle": "round,pad=0.3", "fc": "#FFF0F3", "ec": "#D90429", "lw": 1},
        )

        # Highlight Leu769 -- Leu183
        anchor_dist = _get_pair_dist(769, 183)
        anchor_text = (
            f"Leu769 -- Leu183\n({anchor_dist:.2f} Å Hydrophobic Anchor)"
            if anchor_dist is not None
            else "Leu769 -- Leu183\n(Hydrophobic Anchor)"
        )
        ax1.annotate(
            anchor_text,
            xy=(183, 769),
            xytext=(195, 710),
            arrowprops={"facecolor": "#2B9348", "shrink": 0.08, "width": 1.5, "headwidth": 6},
            fontsize=8.5,
            fontweight="bold",
            bbox={"boxstyle": "round,pad=0.3", "fc": "#E8F5E9", "ec": "#2B9348", "lw": 1},
        )

        # Panel B plotting ...
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        save_figure_safely(fig, output_path)
    finally:
        plt.close(fig)

    return output_path
```

### Step 3: Remediate Figure 3 (`generate_figure3_tpp1_conformational_shifts`)
1. Wrap the entire plotting body in `try: ... finally: plt.close(fig)`.
2. Dynamically extract target coordinates and annotation text from `df_tpp1`:
```python
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(13, 8), sharex=True, gridspec_kw={"height_ratios": [2.2, 1]})
    try:
        # Panel A plotting ...
        disp_map = (
            df_tpp1.set_index("auth_seq_id")["ca_distance_7tre_vs_2i46_angstrom"].to_dict()
            if "auth_seq_id" in df_tpp1.columns and "ca_distance_7tre_vs_2i46_angstrom" in df_tpp1.columns
            else {}
        )

        # Lys170 (Peak)
        d_170 = disp_map.get(170)
        if d_170 is not None:
            ax1.annotate(
                f"Lys170\n(Peak: {d_170:.2f} Å)",
                xy=(170, d_170),
                xytext=(160, max(d_170 + 1.5, 5.8)),
                arrowprops={"facecolor": "#D90429", "shrink": 0.08, "width": 1.5, "headwidth": 6},
                fontsize=8.5,
                fontweight="bold",
                bbox={"boxstyle": "round,pad=0.3", "fc": "#FFF0F3", "ec": "#D90429", "lw": 1},
            )

        # Glu169
        d_169 = disp_map.get(169)
        if d_169 is not None:
            ax1.annotate(
                f"Glu169\n({d_169:.2f} Å)",
                xy=(169, d_169),
                xytext=(145, max(d_169 + 0.7, 4.5)),
                arrowprops={"facecolor": "#D90429", "shrink": 0.08, "width": 1.2, "headwidth": 5},
                fontsize=8,
                bbox={"boxstyle": "round,pad=0.2", "fc": "#FFF0F3", "ec": "#D90429", "lw": 0.8},
            )

        # Glu168
        d_168 = disp_map.get(168)
        if d_168 is not None:
            ax1.annotate(
                f"Glu168\n({d_168:.2f} Å)",
                xy=(168, d_168),
                xytext=(145, max(d_168 + 0.9, 2.2)),
                arrowprops={"facecolor": "#D90429", "shrink": 0.08, "width": 1.2, "headwidth": 5},
                fontsize=8,
                bbox={"boxstyle": "round,pad=0.2", "fc": "#FFF0F3", "ec": "#D90429", "lw": 0.8},
            )

        # Leu183 (Invariant Anchor)
        d_183 = disp_map.get(183)
        if d_183 is not None:
            ax1.annotate(
                f"Leu183 Rigid Anchor\n({d_183:.2f} Å Displacement)",
                xy=(183, d_183),
                xytext=(195, 2.8),
                arrowprops={"facecolor": "#2B9348", "shrink": 0.08, "width": 1.5, "headwidth": 6},
                fontsize=8.5,
                fontweight="bold",
                bbox={"boxstyle": "round,pad=0.3", "fc": "#E8F5E9", "ec": "#2B9348", "lw": 1},
            )

        # Leu212
        d_212 = disp_map.get(212)
        if d_212 is not None:
            ax1.annotate(
                f"Leu212\n({d_212:.2f} Å)",
                xy=(212, d_212),
                xytext=(218, 2.2),
                arrowprops={"facecolor": "#2B9348", "shrink": 0.08, "width": 1.2, "headwidth": 5},
                fontsize=8,
                bbox={"boxstyle": "round,pad=0.2", "fc": "#E8F5E9", "ec": "#2B9348", "lw": 0.8},
            )

        # Dynamic x-axis range
        if len(resnums) > 0:
            ax2.set_xlim(resnums.min() - 2, resnums.max() + 2)

        # Panel B plotting ...
        plt.tight_layout(rect=[0, 0, 1, 0.96])
        save_figure_safely(fig, output_path)
    finally:
        plt.close(fig)

    return output_path
```

### Step 4: Apply `try...finally: plt.close(fig)` & `save_figure_safely` to Figures 2 and 4
In `generate_figure2_cross_state_rmsd`:
```python
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 6))
    try:
        # Plotting ...
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        save_figure_safely(fig, output_path)
    finally:
        plt.close(fig)
```
In `generate_figure4_cutoff_robustness`:
```python
    fig = plt.figure(figsize=(14, 6))
    try:
        # Plotting ...
        plt.tight_layout(rect=[0, 0, 1, 0.95])
        save_figure_safely(fig, output_path)
    finally:
        plt.close(fig)
```

### Step 5: Post-Generation Assertion in `generate_all_figures`
```python
def generate_all_figures(
    tables_dir: Path | str,
    figures_dir: Path | str,
) -> dict[str, Path]:
    tables_dir = Path(tables_dir)
    figures_dir = Path(figures_dir)
    figures_dir.mkdir(parents=True, exist_ok=True)

    fig1_path = figures_dir / "figure1_interface_contacts.png"
    fig2_path = figures_dir / "figure2_cross_state_rmsd.png"
    fig3_path = figures_dir / "figure3_tpp1_conformational_shifts.png"
    fig4_path = figures_dir / "figure4_cutoff_robustness.png"

    generate_figure1_interface_contacts(tables_dir, fig1_path)
    generate_figure2_cross_state_rmsd(tables_dir, fig2_path)
    generate_figure3_tpp1_conformational_shifts(tables_dir, fig3_path)
    generate_figure4_cutoff_robustness(tables_dir, fig4_path)

    figures = {
        "figure1_interface_contacts.png": fig1_path,
        "figure2_cross_state_rmsd.png": fig2_path,
        "figure3_tpp1_conformational_shifts.png": fig3_path,
        "figure4_cutoff_robustness.png": fig4_path,
    }

    # Verify all generated artifacts are present and valid (> 1024 bytes)
    for fname, fpath in figures.items():
        if not fpath.exists():
            raise FileNotFoundError(f"Expected publication figure {fname} was not found at {fpath}")
        size = fpath.stat().st_size
        if size <= 1024:
            raise IOError(f"Figure {fname} is corrupted or empty ({size} bytes, expected > 1024 bytes)")

    return figures
```

---

## 5. Verification Method

To verify these fixes after worker implementation:
1. **Verify Static Code & Patterns**:
   - Check that `plt.close(fig)` is called in a `finally` block in each of the 4 figure functions:
     ```bash
     grep -n "plt.close(fig)" src/telotrade/figures.py
     ```
     Should return 4 matches.
   - Check that `save_figure_safely` is used across all 4 functions:
     ```bash
     grep -n "save_figure_safely" src/telotrade/figures.py
     ```
     Should return definition + 4 calls.
   - Check that `xy=(170, 4.293)` and `xy=(169, 3.809)` are completely gone from Figure 3:
     ```bash
     grep -n "4.293" src/telotrade/figures.py
     ```
     Should return 0 matches.
2. **Linter**:
   ```bash
   ruff check src/ tests/
   ```
   Must pass with 0 errors.
3. **Pipeline Execution & File Sizes**:
   ```bash
   python run_pipeline.py
   ```
   Verify that all 4 figures exist in `results/figures/` and each file size exceeds 50 KB.
4. **Automated Test Suite**:
   ```bash
   pytest tests/ -v
   ```
   All tests, including `test_t4_figures_directory_and_files`, must pass cleanly (100% pass rate).
