import numpy as np

# Standard Bondi Van der Waals radii (in Angstroms)
VDW_RADII = {
    'C': 1.70,
    'N': 1.55,
    'O': 1.52,
    'S': 1.80,
    'P': 1.80,
    'H': 1.20,
    'DEFAULT': 1.70
}

def fibonacci_sphere(samples=100):
    points = []
    phi = np.pi * (3.0 - np.sqrt(5.0))  # golden angle in radians
    for i in range(samples):
        y = 1.0 - (i / float(samples - 1)) * 2.0  # y goes from 1 to -1
        radius = np.sqrt(1.0 - y * y)  # radius at y
        theta = phi * i  # golden angle increment
        x = np.cos(theta) * radius
        z = np.sin(theta) * radius
        points.append([x, y, z])
    return np.array(points, dtype=np.float32)

def compute_sasa(coords, elements, probe_radius=1.4, n_sphere_points=100):
    n_atoms = len(coords)
    sphere_pts = fibonacci_sphere(n_sphere_points)
    
    atom_vdw = np.array([VDW_RADII.get(el, VDW_RADII['DEFAULT']) for el in elements], dtype=np.float32)
    atom_radii = atom_vdw + probe_radius
    
    # Grid hashing for neighbor lookup
    max_radius = np.max(atom_radii)
    cell_size = 2.0 * max_radius
    grid = {}
    for idx, pt in enumerate(coords):
        cell = (int(pt[0] // cell_size), int(pt[1] // cell_size), int(pt[2] // cell_size))
        if cell not in grid:
            grid[cell] = []
        grid[cell].append(idx)
        
    sasa_per_atom = np.zeros(n_atoms, dtype=np.float32)
    
    for i in range(n_atoms):
        pt_i = coords[i]
        r_i = atom_radii[i]
        cx, cy, cz = int(pt_i[0] // cell_size), int(pt_i[1] // cell_size), int(pt_i[2] // cell_size)
        
        # Find potential neighbors
        neighbors = []
        for dx in (-1, 0, 1):
            for dy in (-1, 0, 1):
                for dz in (-1, 0, 1):
                    neighbor_cell = (cx + dx, cy + dy, cz + dz)
                    if neighbor_cell in grid:
                        for j in grid[neighbor_cell]:
                            if i != j:
                                d_sq = np.sum((pt_i - coords[j])**2)
                                if d_sq < (r_i + atom_radii[j])**2:
                                    neighbors.append(j)
        
        if not neighbors:
            # Completely exposed
            sasa_per_atom[i] = 4.0 * np.pi * (r_i**2)
            continue
            
        test_pts = pt_i + sphere_pts * r_i # shape (n_sphere_points, 3)
        accessible = np.ones(n_sphere_points, dtype=bool)
        
        nbr_coords = coords[neighbors]
        nbr_radii_sq = (atom_radii[neighbors])**2
        
        for k in range(len(neighbors)):
            diff = test_pts - nbr_coords[k]
            d2 = np.sum(diff**2, axis=1)
            accessible[d2 < nbr_radii_sq[k]] = False
            if not np.any(accessible):
                break
                
        frac_accessible = np.sum(accessible) / float(n_sphere_points)
        sasa_per_atom[i] = 4.0 * np.pi * (r_i**2) * frac_accessible
        
    return sasa_per_atom

print("Testing Fibonacci sphere Shrake-Rupley generator...")
pts = fibonacci_sphere(100)
print("Points shape:", pts.shape)
print("Norm check (first 5):", np.linalg.norm(pts[:5], axis=1))
