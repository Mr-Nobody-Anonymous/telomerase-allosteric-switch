import os
import json
from parse_2i46_seq import parse_simple_cif_loop

pdb_ids = ['7TRE', '7TRD', '7TRF', '7TRC', '2I46']

for pid in pdb_ids:
    cif_path = f'cifs/{pid}.cif'
    with open(cif_path, 'r', encoding='utf-8') as f:
        lines = f.readlines()
    
    print('======================================================================')
    print(f'PDB ID: {pid}')
    
    # Check exptl
    for l in lines:
        if l.startswith('_exptl.method'):
            print('  Exptl method:', l.strip())
        if l.startswith('_refine.ls_d_res_high'):
            print('  _refine.ls_d_res_high:', l.strip())
        if l.startswith('_reflns.d_resolution_high'):
            print('  _reflns.d_resolution_high:', l.strip())
        if l.startswith('_em_3d_reconstruction.resolution'):
            print('  EM 3D reconstruction resolution:', l.strip())
        if l.startswith('_em_3d_reconstruction.method'):
            print('  EM 3D reconstruction method:', l.strip())
        if l.startswith('_em_3d_reconstruction.symmetry_type'):
            print('  EM symmetry:', l.strip())
        if l.startswith('_em_imaging.microscope_model'):
            print('  Microscope:', l.strip())
        if l.startswith('_em_imaging.accelerating_voltage'):
            print('  Voltage:', l.strip())
        if l.startswith('_em_imaging.detector_mode'):
            print('  Detector mode:', l.strip())
        if l.startswith('_em_detector.detector'):
            print('  Detector:', l.strip())
        if l.startswith('_refine.ls_R_factor_R_work'):
            print('  R_work:', l.strip())
        if l.startswith('_refine.ls_R_factor_R_free'):
            print('  R_free:', l.strip())
        if l.startswith('_refine.B_iso_mean') or l.startswith('_refine.B_iso_max'):
            print('  B_factor:', l.strip())
        if l.startswith('_em_3d_reconstruction.fsc_cutoff'):
            print('  FSC cutoff:', l.strip())

    # Check validation report details from RCSB JSON if available
    # Also check EMDB ID
    for l in lines:
        if 'EMD-' in l:
            print('  EMDB reference in CIF:', l.strip()[:80])
            break
