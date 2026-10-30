import gemmi
from openmm.app import PDBFile
from pdbfixer import PDBFixer

data_dir = "C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data"

print("1. Preparing Apo WT fixed structure...")
fixer_apo = PDBFixer(filename=f"{data_dir}/junction_apo_continuous.pdb")
fixer_apo.findMissingResidues()
fixer_apo.findMissingAtoms()
fixer_apo.addMissingAtoms()
fixer_apo.addMissingHydrogens(7.5)
with open(f"{data_dir}/junction_apo_fixed.pdb", 'w') as f_out:
    PDBFile.writeFile(fixer_apo.topology, fixer_apo.positions, f_out)
print("   Apo WT fixed:", fixer_apo.topology.getNumAtoms(), "atoms.")

print("2. Generating K757A mutant from Bound WT...")
# Load 7TRE bound fixed PDB with gemmi to mutate Lys757
st_bound = gemmi.read_structure(f"{data_dir}/junction_bound_fixed.pdb")
# Find Lys757 in chain B (TRAP)
for r in st_bound[0]["B"]:
    if r.seqid.num == 757:
        r.name = "ALA"
        # Delete sidechain atoms beyond CB
        atoms_to_del = [a.name for a in r if a.name in ["CG", "CD", "CE", "NZ", "HZ1", "HZ2", "HZ3", "1HZ", "2HZ", "3HZ", "HD2", "HD3", "HE2", "HE3", "HG2", "HG3"]]
        for aname in atoms_to_del:
            del r[aname]
st_bound.write_pdb(f"{data_dir}/junction_bound_k757a_raw.pdb")

# Run pdbfixer to add correct alanine hydrogens
fixer_k757a = PDBFixer(filename=f"{data_dir}/junction_bound_k757a_raw.pdb")
fixer_k757a.findMissingResidues()
fixer_k757a.findMissingAtoms()
fixer_k757a.addMissingAtoms()
fixer_k757a.addMissingHydrogens(7.5)
with open(f"{data_dir}/junction_bound_k757a_fixed.pdb", 'w') as f_out:
    PDBFile.writeFile(fixer_k757a.topology, fixer_k757a.positions, f_out)
print("   Bound K757A fixed:", fixer_k757a.topology.getNumAtoms(), "atoms.")

print("All MD input systems prepared successfully!")
