import numpy as np
from openmm import LangevinMiddleIntegrator, Platform
from openmm.app import ForceField, NoCutoff, PDBFile, Simulation
from openmm.unit import femtoseconds, kelvin, nanometers, picoseconds
from scipy import stats

data_dir = "C:/Users/hp/.gemini/antigravity/scratch/telomerase-research/data"
ff = ForceField('amber14-all.xml', 'implicit/obc2.xml')
platform = Platform.getPlatformByName('CPU')

def analyze_replicate_with_burnin(pdb_path, seed, n_steps=1500, report_interval=20, burnin_steps=300):
    pdb = PDBFile(pdb_path)
    system = ff.createSystem(pdb.topology, nonbondedMethod=NoCutoff)
    integrator = LangevinMiddleIntegrator(300*kelvin, 1/picoseconds, 2*femtoseconds)
    integrator.setRandomNumberSeed(seed)
    
    sim = Simulation(pdb.topology, system, integrator, platform)
    sim.context.setPositions(pdb.positions)
    sim.minimizeEnergy(maxIterations=150)
    
    # Map CA atoms and target atoms
    top = pdb.topology
    atom_y_oh = None
    atom_q_oe1 = None
    ca_atoms = []
    
    for r in top.residues():
        c_id = r.chain.id
        for a in r.atoms():
            if a.name == 'CA':
                ca_atoms.append(a.index)
            if c_id == 'A' and r.name == 'TYR' and a.name == 'OH':
                atom_y_oh = a.index
            if c_id == 'B' and r.name == 'GLN' and a.name == 'OE1':
                atom_q_oe1 = a.index
                
    init_ca = np.array(pdb.positions.value_in_unit(nanometers))[ca_atoms] * 10.0 # initial CA coords in A
    
    rmsd_traj = []
    dist_traj = []
    steps_list = []
    
    for step in range(0, n_steps, report_interval):
        sim.step(report_interval)
        state = sim.context.getState(getPositions=True)
        pos = state.getPositions(asNumpy=True).value_in_unit(nanometers) * 10.0
        
        # Compute CA RMSD relative to starting structure
        curr_ca = pos[ca_atoms]
        # Centering for RMSD
        c_curr = curr_ca - np.mean(curr_ca, axis=0)
        c_init = init_ca - np.mean(init_ca, axis=0)
        rmsd = np.sqrt(np.mean(np.sum((c_curr - c_init)**2, axis=-1)))
        rmsd_traj.append(rmsd)
        
        # Interdomain distance
        if atom_y_oh is not None and atom_q_oe1 is not None:
            dist = np.linalg.norm(pos[atom_y_oh] - pos[atom_q_oe1])
            dist_traj.append(dist)
            
        steps_list.append(step + report_interval)
        
    rmsd_traj = np.array(rmsd_traj)
    dist_traj = np.array(dist_traj)
    steps_list = np.array(steps_list)
    
    # Separate burn-in from production phase
    prod_mask = steps_list > burnin_steps
    prod_dist = dist_traj[prod_mask]
    prod_rmsd = rmsd_traj[prod_mask]
    
    mean_prod_dist = np.mean(prod_dist)
    std_prod_dist = np.std(prod_dist)
    hbond_occ_prod = np.mean(prod_dist <= 3.2) * 100.0
    mean_rmsd_prod = np.mean(prod_rmsd)
    
    return {
        "steps": steps_list,
        "rmsd": rmsd_traj,
        "dist": dist_traj,
        "mean_prod_dist": mean_prod_dist,
        "std_prod_dist": std_prod_dist,
        "hbond_occ_prod": hbond_occ_prod,
        "mean_rmsd_prod": mean_rmsd_prod
    }

print("Running rigorous equilibration audit across N=3 replicates with 300-step burn-in...")
seeds_bound = [101, 102, 103]
seeds_apo = [201, 202, 203]

bound_data = []
for s in seeds_bound:
    res = analyze_replicate_with_burnin(f"{data_dir}/junction_bound_fixed.pdb", seed=s)
    bound_data.append(res)
    print(f"  Bound WT (Seed {s}): RMSD plateau = {res['mean_rmsd_prod']:.2f} A | Post-burnin Dist = {res['mean_prod_dist']:.3f} +/- {res['std_prod_dist']:.3f} A | H-bond Occ = {res['hbond_occ_prod']:.1f}%")

apo_data = []
for s in seeds_apo:
    res = analyze_replicate_with_burnin(f"{data_dir}/junction_apo_fixed.pdb", seed=s)
    apo_data.append(res)
    print(f"  Apo WT (Seed {s}): RMSD plateau = {res['mean_rmsd_prod']:.2f} A | Post-burnin Dist = {res['mean_prod_dist']:.3f} +/- {res['std_prod_dist']:.3f} A | H-bond Occ = {res['hbond_occ_prod']:.1f}%")

bound_means = [d['mean_prod_dist'] for d in bound_data]
apo_means = [d['mean_prod_dist'] for d in apo_data]
bound_occs = [d['hbond_occ_prod'] for d in bound_data]
apo_occs = [d['hbond_occ_prod'] for d in apo_data]

t_d, p_d = stats.ttest_ind(bound_means, apo_means)
t_o, p_o = stats.ttest_ind(bound_occs, apo_occs)

print("\n--- EQUILIBRATED PRODUCTION STATISTICAL COMPARISON (N=3 Replicates as Statistical Units) ---")
print(f"Equilibrated Distance: Bound = {np.mean(bound_means):.3f} +/- {np.std(bound_means):.3f} A vs Apo = {np.mean(apo_means):.3f} +/- {np.std(apo_means):.3f} A")
print(f"Two-sample t-test on independent replicas: t = {t_d:.3f}, p-value = {p_d:.4e}")
print(f"Equilibrated H-Bond Occupancy: Bound = {np.mean(bound_occs):.1f}% +/- {np.std(bound_occs):.1f}% vs Apo = {np.mean(apo_occs):.1f}% +/- {np.std(apo_occs):.1f}%")
print(f"Two-sample t-test on independent replicas: t = {t_o:.3f}, p-value = {p_o:.4e}")
