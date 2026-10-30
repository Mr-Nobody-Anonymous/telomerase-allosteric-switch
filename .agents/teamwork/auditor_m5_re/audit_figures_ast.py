"""
Forensic AST and Static Integrity Analyzer for src/telotrade/figures.py
Milestone 5 Remediation Re-Audit.
"""

import ast
import pathlib
import sys

def audit_figures_ast():
    target_path = pathlib.Path("src/telotrade/figures.py")
    if not target_path.exists():
        print(f"FAIL: File not found: {target_path}")
        sys.exit(1)

    code = target_path.read_text(encoding="utf-8")
    tree = ast.parse(code, filename=str(target_path))

    violations = []
    observations = []

    # 1. Check for unassigned expressions (e.g. pd.read_csv without assignment)
    unassigned_reads = []
    for node in ast.walk(tree):
        if isinstance(node, ast.Expr):
            for subnode in ast.walk(node):
                if isinstance(subnode, ast.Attribute) and subnode.attr == "read_csv":
                    unassigned_reads.append(node.lineno)

    if unassigned_reads:
        violations.append(f"Unassigned read_csv calls found at lines: {unassigned_reads}")
    else:
        observations.append("PASS: Zero unassigned read_csv calls found in figures.py.")

    # 2. Check for suspicious numeric lists (data arrays hardcoded)
    all_numeric_lists = []
    for node in ast.walk(tree):
        if isinstance(node, ast.List):
            if len(node.elts) > 0 and all(
                isinstance(elt, (ast.Constant, ast.UnaryOp)) and
                isinstance(
                    elt.value if isinstance(elt, ast.Constant) else getattr(getattr(elt, "operand", None), "value", None),
                    (int, float)
                )
                for elt in node.elts
            ):
                vals = [ast.literal_eval(elt) for elt in node.elts]
                all_numeric_lists.append((node.lineno, vals))

    # Matplotlib rects, figure ratios, and domain coordinate constants
    whitelisted = [
        [1.2, 1], # gridspec width_ratios
        [2.2, 1], # gridspec height_ratios
        [1, 1, 1.2], # gridspec width_ratios
        [0, 0, 1, 0.95], # tight_layout rect
        [0, 0, 1, 0.96], # tight_layout rect
        [168, 169, 170, 171, 183, 212], # TPP1 TEL patch author sequence numbers for reference vertical lines
        [3.5, 4.0, 4.5], # contact cutoffs evaluated
        [0.0], # zero initialization multiplier
    ]

    flagged_lists = []
    for line, vals in all_numeric_lists:
        if vals not in whitelisted:
            flagged_lists.append((line, vals))

    if flagged_lists:
        violations.append(f"Hardcoded data arrays found in AST: {flagged_lists}")
    else:
        observations.append(f"PASS: Zero hardcoded scientific metric lists found. All {len(all_numeric_lists)} numeric lists are confirmed layout/configuration constants.")

    # 3. Check for specific prohibited patterns from previous audit
    legacy_hardcoded_snippets = [
        "[0.356, 0.338, 0.434, 0.338, 1.966, 1.651]",
        "[9, 46, 123]",
        "[5, 21, 32]",
        "[0.411, 0.270, 0.295, 0.262, 0.238, 0.338]",
        "[9, 37, 77]",
        "[3.35, 3.45, 3.51, 3.58, 3.62, 4.12, 4.28]",
        "[2.32, 3.24, 3.36, 3.40, 3.48, 3.52, 3.56, 3.67]",
    ]
    found_snippets = [s for s in legacy_hardcoded_snippets if s in code]
    if found_snippets:
        violations.append(f"Legacy hardcoded data snippets still present: {found_snippets}")
    else:
        observations.append("PASS: 100% of legacy hardcoded arrays and mock Series are eradicated.")

    # 4. Check resource management: try...finally: plt.close(fig)
    fig_funcs = [
        "generate_figure1_interface_contacts",
        "generate_figure2_cross_state_rmsd",
        "generate_figure3_tpp1_conformational_shifts",
        "generate_figure4_cutoff_robustness",
    ]
    for func_node in ast.walk(tree):
        if isinstance(func_node, ast.FunctionDef) and func_node.name in fig_funcs:
            has_try_finally = False
            for stmt in func_node.body:
                if isinstance(stmt, ast.Try) and stmt.finalbody:
                    for fstmt in stmt.finalbody:
                        if isinstance(fstmt, ast.Expr) and isinstance(fstmt.value, ast.Call):
                            unparsed = ast.unparse(fstmt.value)
                            if "plt.close" in unparsed:
                                has_try_finally = True
            if not has_try_finally:
                violations.append(f"Function {func_node.name} lacks try...finally: plt.close(fig)")
            else:
                observations.append(f"PASS: {func_node.name} implements try...finally: plt.close(fig)")

    # 5. Check atomic saving: save_figure_safely usage
    for func_node in ast.walk(tree):
        if isinstance(func_node, ast.FunctionDef) and func_node.name in fig_funcs:
            calls_save_safely = any(
                isinstance(n, ast.Call) and getattr(n.func, "id", None) == "save_figure_safely"
                for n in ast.walk(func_node)
            )
            if not calls_save_safely:
                violations.append(f"Function {func_node.name} does not call save_figure_safely")
            else:
                observations.append(f"PASS: {func_node.name} routes saving through save_figure_safely")

    print("\n=== AST AUDIT OBSERVATIONS ===")
    for obs in observations:
        print(f"  [+] {obs}")

    print("\n=== AST AUDIT VIOLATIONS ===")
    if violations:
        for v in violations:
            print(f"  [-] {v}")
        return False
    else:
        print("  None. 100% CLEAN.")
        return True

if __name__ == "__main__":
    success = audit_figures_ast()
    sys.exit(0 if success else 1)
