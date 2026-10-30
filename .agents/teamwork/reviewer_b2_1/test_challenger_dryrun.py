import pytest
import sys
from pathlib import Path

# Run pytest specifically on tests/test_adversarial_b2_challenger.py, but dynamically patch Bio.PDB.ShrakeRupley
import Bio.PDB
import types

mod = types.ModuleType("Bio.PDB.ShrakeRupley")
mod.ShrakeRupley = Bio.PDB.ShrakeRupley
sys.modules["Bio.PDB.ShrakeRupley"] = mod

code = pytest.main(["tests/test_adversarial_b2_challenger.py", "-v"])
print("Exit code:", code)
