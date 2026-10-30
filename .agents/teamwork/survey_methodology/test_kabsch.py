import numpy as np

def kabsch_fit(p, q):
    """
    Superimposes coordinates P onto Q using the Kabsch algorithm (SVD).
    p, q: numpy arrays of shape (N, 3)
    Returns:
      rmsd: float
      r: rotation matrix (3, 3)
      t: translation vector (3,)
      p_aligned: transformed p coordinates (N, 3)
    """
    assert p.shape == q.shape, "Shapes must match"
    n = p.shape[0]
    centroid_p = np.mean(p, axis=0)
    centroid_q = np.mean(q, axis=0)
    
    p_c = p - centroid_p
    q_c = q - centroid_q
    
    h = np.dot(p_c.T, q_c)
    u, s, vt = np.linalg.svd(h)
    
    # Check reflection
    d = np.linalg.det(np.dot(vt.T, u.T))
    e = np.eye(3)
    if d < 0:
        e[2, 2] = -1.0
        
    r = np.dot(vt.T, np.dot(e, u.T))
    t = centroid_q - np.dot(r, centroid_p)
    
    p_aligned = np.dot(p, r.T) + t
    rmsd = np.sqrt(np.mean(np.sum((p_aligned - q)**2, axis=1)))
    return rmsd, r, t, p_aligned

# Simple test
np.random.seed(42)
pts = np.random.randn(50, 3)
# Rotate and translate
rot = np.array([[0, -1, 0], [1, 0, 0], [0, 0, 1]], dtype=float)
pts_q = np.dot(pts, rot.T) + np.array([10.0, 5.0, -2.0])

rmsd, r, t, p_aligned = kabsch_fit(pts, pts_q)
print(f"Kabsch test RMSD on exact rigid transform: {rmsd:.6e}")
print("Recovered rotation close to true?", np.allclose(r, rot))
