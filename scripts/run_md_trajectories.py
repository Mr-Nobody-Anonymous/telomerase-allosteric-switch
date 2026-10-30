import numpy as np
from openmm import LangevinMiddleIntegrator, Platform
from openmm.app import ForceField, NoCutoff, PDBFile, Simulation
from openmm.unit import femtoseconds, kelvin, nanometers, picoseconds
from scipy import stats

data_dir = "C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data"
ff = ForceField('amber14-all.xml', 'implicit/obc2.xml')
platform = Platform.getPlatformByName('CPU')

def run_simulation_replicate(pdb_path, seed, n_steps=1000, report_interval=20):
    pdb = PDBFile(pdb_path)
    system = ff.createSystem(pdb.topology, nonbondedMethod=NoCutoff)
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picoseconds, 2*femtoseconds)
    integrator.setRandomNumberSeed(seed)
    
    sim = Simulation(pdb.topology, system, integrator, platform)
    sim.context.setPositions(pdb.positions)
    sim.minimizeEnergy(maxIterations=100)
    
    # Identify atom indices
    # Chain A residue 11 is Tyr176
    # Chain B residue 44 is Gln794
    # Chain B residue 7 is Lys757
    top = pdb.topology
    res_list = list(top.residues())
    
    atom_y_oh = None
    atom_q_oe1 = None
    ca_atoms = []
    ca_res_map = []
    
    for r in res_list:
        c_id = r.chain.id
        for a in r.atoms():
            if a.name == 'CA':
                ca_atoms.append(a.index)
                ca_res_map.append((c_id, r.name, r.id))
            if c_id == 'A' and r.name == 'TYR' and a.name == 'OH':
                atom_y_oh = a.index
            if c_id == 'B' and r.name == 'GLN' and a.name == 'OE1':
                atom_q_oe1 = a.index
                
    distances = []
    ca_coords_traj = []
    
    for step in range(0, n_steps, report_interval):
        sim.step(report_interval)
        state = sim.context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(nanometers) * 10.0 # convert nm to A
        
        # Distance Tyr176:OH to Gln794:OE1
        if atom_y_oh is not None and atom_q_oe1 is not None:
            p1 = pos[atom_y_oh]
            p2 = pos[atom_q_oe1]
            dist = np.linalg.norm(p1 - p2)
            distances.append(dist)
            
        ca_coords_traj.append(pos[ca_atoms])
        
    ca_coords_traj = np.array(ca_coords_traj) # shape: (frames, n_ca, 3)
    # Calculate per-residue RMSF
    mean_coords = np.mean(ca_coords_traj, axis=0)
    rmsf = np.sqrt(np.mean(np.sum((ca_coords_traj - mean_coords)**2, axis=-1), axis=0))
    
    hbond_occupancy = np.mean(np.array(distances) <= 3.2) * 100.0
    mean_dist = np.mean(distances)
    std_dist = np.std(distances)
    
    return {
        "distances": distances,
        "mean_dist": mean_dist,
        "std_dist": std_dist,
        "hbond_occ": hbond_occupancy,
        "rmsf": rmsf,
        "ca_map": ca_res_map
    }

print("=================================================================")
print("EXECUTING OPENMM ALL-ATOM MOLECULAR DYNAMICS SIMULATION SUITE")
print("=================================================================")

seeds_bound = [101, 102, 103]
seeds_apo = [201, 202, 203]

print("\n--- 1. Simulating TPP1-Bound WT (3 Replicates) ---")
bound_results = []
for s in seeds_bound:
    print(f"  Running Bound WT Replicate Seed {s}...")
    res = run_simulation_replicate(f"{data_dir}/junction_bound_fixed.pdb", seed=s)
    bound_results.append(res)
    print(f"    Seed {s}: Mean H-bond dist = {res['mean_dist']:.3f} +/- {res['std_dist']:.3f} A | H-bond occupancy = {res['hbond_occ']:.1f}%")

print("\n--- 2. Simulating Apo WT (3 Replicates) ---")
apo_results = []
for s in seeds_apo:
    print(f"  Running Apo WT Replicate Seed {s}...")
    res = run_simulation_replicate(f"{data_dir}/junction_apo_fixed.pdb", seed=s)
    apo_results.append(res)
    print(f"    Seed {s}: Mean H-bond dist = {res['mean_dist']:.3f} +/- {res['std_dist']:.3f} A | H-bond occupancy = {res['hbond_occ']:.1f}%")

# Aggregate statistics
bound_dists = [r['mean_dist'] for r in bound_results]
apo_dists = [r['mean_dist'] for r in apo_results]
bound_occs = [r['hbond_occ'] for r in bound_results]
apo_occs = [r['hbond_occ'] for r in apo_results]

t_stat_dist, p_val_dist = stats.ttest_ind(bound_dists, apo_dists)
t_stat_occ, p_val_occ = stats.ttest_ind(bound_occs, apo_occs)

print("\n=================================================================")
print("STATISTICAL SUMMARY ACROSS INDEPENDENT MD REPLICATES (N=3)")
print("=================================================================")
print("Tyr176:OH <-> Gln794:OE1 Inter-Domain Distance:")
print(f"  Bound WT: {np.mean(bound_dists):.3f} +/- {np.std(bound_dists):.3f} A")
print(f"  Apo WT:   {np.mean(apo_dists):.3f} +/- {np.std(apo_dists):.3f} A")
print(f"  Two-sample t-test: t = {t_stat_dist:.3f}, p-value = {p_val_dist:.4e}")

print("\nHydrogen-Bond Occupancy (<= 3.2 A):")
print(f"  Bound WT: {np.mean(bound_occs):.1f}% +/- {np.std(bound_occs):.1f}%")
print(f"  Apo WT:   {np.mean(apo_occs):.1f}% +/- {np.std(apo_occs):.1f}%")
print(f"  Two-sample t-test: t = {t_stat_occ:.3f}, p-value = {p_val_occ:.4e}")

# Save results for markdown artifact
with open(f"{data_dir}/md_replicate_summary.txt", "w") as f:
    f.write(f"Bound_Mean_Dist={np.mean(bound_dists):.4f}\n")
    f.write(f"Bound_Std_Dist={np.std(bound_dists):.4f}\n")
    f.write(f"Apo_Mean_Dist={np.mean(apo_dists):.4f}\n")
    f.write(f"Apo_Std_Dist={np.std(apo_dists):.4f}\n")
    f.write(f"P_Val_Dist={p_val_dist:.6e}\n")
    f.write(f"Bound_Mean_Occ={np.mean(bound_occs):.2f}\n")
    f.write(f"Apo_Mean_Occ={np.mean(apo_occs):.2f}\n")
    f.write(f"P_Val_Occ={p_val_occ:.6e}\n")
print("Saved summary metrics to data/md_replicate_summary.txt")
