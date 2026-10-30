import Bio.PDB
import copy

parser = Bio.PDB.MMCIFParser(QUIET=True)
structure = parser.get_structure('7TRE', r'C:\Users\hp\.gemini\antigravity\scratch\telomerase-research\.agents\teamwork\survey_pdb\cifs\7TRE.cif')
model = structure[0]

# Chain A is TERT, Chain C is TPP1
chain_A = model['A']
chain_C = model['C']

sr = Bio.PDB.ShrakeRupley(probe_radius=1.4, n_points=960)

# Compute complex SASA (with just A and C, or full complex)
# To get BSA of the A:C interface specifically:
# SASA_A: isolate chain A alone
# SASA_C: isolate chain C alone
# SASA_AC: complex of chain A + chain C alone (excluding RNA and DNA)

struct_AC = Bio.PDB.Structure.Structure('AC')
model_AC = Bio.PDB.Model.Model(0)
model_AC.add(copy.deepcopy(chain_A))
model_AC.add(copy.deepcopy(chain_C))
struct_AC.add(model_AC)

struct_A = Bio.PDB.Structure.Structure('A')
model_A = Bio.PDB.Model.Model(0)
model_A.add(copy.deepcopy(chain_A))
struct_A.add(model_A)

struct_C = Bio.PDB.Structure.Structure('C')
model_C = Bio.PDB.Model.Model(0)
model_C.add(copy.deepcopy(chain_C))
struct_C.add(model_C)

sr.compute(struct_A, level='A')
sasa_A = sum(atom.sasa for atom in struct_A.get_atoms())

sr.compute(struct_C, level='A')
sasa_C = sum(atom.sasa for atom in struct_C.get_atoms())

sr.compute(struct_AC, level='A')
sasa_AC = sum(atom.sasa for atom in struct_AC.get_atoms())

delta_sasa = sasa_A + sasa_C - sasa_AC
bsa = delta_sasa / 2.0

print(f"SASA TERT (Chain A alone): {sasa_A:.2f} A^2")
print(f"SASA TPP1 (Chain C alone): {sasa_C:.2f} A^2")
print(f"SASA Complex (TERT + TPP1): {sasa_AC:.2f} A^2")
print(f"Delta SASA (Total Buried Area): {delta_sasa:.2f} A^2")
print(f"BSA (Formula: (SASA_A + SASA_C - SASA_AC) / 2): {bsa:.2f} A^2")
