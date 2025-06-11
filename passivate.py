#!/usr/bin/env python3
"""
attach_ligands_cat.py – CAT §2.5 attachment for CsPbBr3 QD + Lecithin

Usage example:
  python attach_ligands_cat.py \
    --qd csrbpbbr3_ordered.xyz \
    --ligand Lecithin.xyz \
    --n_idx 2 --c_idx 109 \
    --offset_out 1.5 \
    --angle_plane 90 \
    --coarse_step_deg 15 \
    --refine_step_deg 1 \
    --warn_tol 1.4 \
    --neigh 3 \
    --out qd_cat_attached.xyz
"""
import argparse, math, numpy as np
from ase import Atoms, io
from scipy.spatial import KDTree

# ——— Geometry helpers ——————————————————————————————————————————————————

def unit(v):
    v = np.asarray(v, float)
    n = np.linalg.norm(v)
    return v/n if n else v

def facet_normal(qd_pos, idx, k):
    kdt = KDTree(qd_pos)
    _, nbrs = kdt.query(qd_pos[idx], k=k+1)
    nbrs = nbrs[1:]
    v1 = qd_pos[nbrs[1]] - qd_pos[nbrs[0]]
    v2 = qd_pos[nbrs[2]] - qd_pos[nbrs[0]]
    return unit(np.cross(v1, v2))

def outward_normal(qd, idx, n):
    com = qd.get_center_of_mass()
    to_anchor = qd[idx].position - com
    return -n if np.dot(to_anchor, n) < 0 else n

def rotation_matrix_from_vectors(a, b):
    a, b = unit(a), unit(b)
    v = np.cross(a, b); s = np.linalg.norm(v)
    if s < 1e-8:
        return np.eye(3) if np.dot(a, b) > 0 else -np.eye(3)
    c = np.dot(a, b)
    vx = np.array([[   0, -v[2],  v[1]],
                   [ v[2],    0, -v[0]],
                   [-v[1], v[0],    0]])
    return np.eye(3) + vx + vx @ vx * ((1 - c) / (s**2))

def make_ligand_template(lig, n_idx, c_idx):
    pos = lig.get_positions().copy()
    vCN = pos[n_idx] - pos[c_idx]
    pos -= pos[n_idx]
    return pos, unit(vCN), lig.get_atomic_numbers()

def rotate_about_axis(coords, axis, ang):
    axis = unit(axis)
    c, s = math.cos(ang), math.sin(ang)
    ux, uy, uz = axis
    R = np.array([
        [c+ux*ux*(1-c),   ux*uy*(1-c)-uz*s, ux*uz*(1-c)+uy*s],
        [uy*ux*(1-c)+uz*s, c+uy*uy*(1-c),   uy*uz*(1-c)-ux*s],
        [uz*ux*(1-c)-uy*s, uz*uy*(1-c)+ux*s, c+uz*uz*(1-c)]
    ])
    return (R @ coords.T).T

def build_heavy_kdtree(qd, placed):
    core_pos = np.array([a.position for a in qd if a.symbol!='Rb' and a.number>1])
    heavy = [core_pos]
    for lig in placed:
        pos = lig.get_positions()
        nums= lig.get_atomic_numbers()
        heavy.append(pos[nums>1])
    return KDTree(np.vstack(heavy))

# ——— CAT §2.5 attachment ——————————————————————————————————————————————————

def place_one_cat(idx, anchor_pos, qd, qd_pos,
                  ref_pos, ref_vec, at_nums,
                  n0, offset_out,
                  core_and_placed,
                  coarse_step_deg, refine_step_deg, warn_tol):
    R_align = rotation_matrix_from_vectors(ref_vec, n0)
    aligned = (R_align @ ref_pos.T).T
    kdtree = build_heavy_kdtree(*core_and_placed)
    heavy_idx = np.where(at_nums>1)[0]

    # coarse scan
    phis = np.deg2rad(np.arange(0,360,coarse_step_deg))
    dmins = []
    for φ in phis:
        posφ = rotate_about_axis(aligned, n0, φ) + anchor_pos + n0*offset_out
        d, _ = kdtree.query(posφ[heavy_idx], k=1)
        dmins.append(d.min())
    i0 = int(np.argmax(dmins))
    best_phi, best_d = phis[i0], dmins[i0]

    # refine ±coarse_step
    for δ in np.deg2rad(np.arange(-coarse_step_deg, coarse_step_deg+refine_step_deg, refine_step_deg)):
        φ = best_phi + δ
        posφ = rotate_about_axis(aligned, n0, φ) + anchor_pos + n0*offset_out
        d, _ = kdtree.query(posφ[heavy_idx], k=1)
        dmin = d.min()
        if dmin>best_d:
            best_d, best_phi = dmin, φ

    final_pos = rotate_about_axis(aligned, n0, best_phi) + anchor_pos + n0*offset_out
    lig = Atoms(numbers=at_nums, positions=final_pos)

    deg = (best_phi*180/math.pi) % 360
    msg = f"    → φₒₚₜ = {deg:.1f}°, min-dist = {best_d:.3f} Å"
    if best_d < warn_tol: msg += f"  ⚠️ (<{warn_tol} Å)"
    print(msg)
    return lig

# ——— Main ——————————————————————————————————————————————————————————

def main():
    p = argparse.ArgumentParser()
    p.add_argument('--qd',             default='csrbpbbr3_ordered.xyz')
    p.add_argument('--ligand',         default='Lecithin.xyz')
    p.add_argument('--n_idx',          type=int, default=2)
    p.add_argument('--c_idx',          type=int, default=109)
    p.add_argument('--offset_out',     type=float, default=1.5)
    p.add_argument('--angle_plane',    type=float, default=None)
    p.add_argument('--coarse_step_deg',type=float, default=15.0)
    p.add_argument('--refine_step_deg',type=float, default=1.0)
    p.add_argument('--warn_tol',       type=float, default=1.4)
    p.add_argument('--neigh',          type=int, default=3)
    p.add_argument('--out',            default='qd_cat_attached.xyz')
    args = p.parse_args()

    if args.angle_plane is not None:
        tilt = abs(args.angle_plane - 90.0)
        print(f"→ angle_plane={args.angle_plane}° ⇒ tilt-from-normal={tilt:.1f}°")
    else:
        tilt = 0.0
        print("→ default tilt-from-normal = 0.0°")

    print("Parameters:")
    for k in ('qd','ligand','n_idx','c_idx','offset_out',
              'angle_plane','coarse_step_deg','refine_step_deg',
              'warn_tol','neigh','out'):
        print(f"  {k:14s}: {getattr(args,k)}")
    print()

    qd     = io.read(args.qd)
    ligand = io.read(args.ligand)
    print(f"Loaded QD ({len(qd)} atoms) and ligand ({len(ligand)} atoms).")

    n_idx, c_idx = args.n_idx-1, args.c_idx-1
    ref_pos, ref_vec, at_nums = make_ligand_template(ligand, n_idx, c_idx)
    qd_pos = qd.get_positions()
    anchors = [i for i,a in enumerate(qd) if a.symbol=='Rb']
    print(f"Found {len(anchors)} Rb anchors.\n")

    placed = []
    for i, idx in enumerate(anchors,1):
        pos = qd_pos[idx]
        print(f"Anchor {i}/{len(anchors)} (atom {idx} @ {pos}):")
        n0 = facet_normal(qd_pos, idx, args.neigh)
        n0 = outward_normal(qd, idx, n0)
        lig_i = place_one_cat(
            idx, pos, qd, qd_pos,
            ref_pos, ref_vec, at_nums,
            n0, args.offset_out,
            (qd, placed),
            args.coarse_step_deg,
            args.refine_step_deg,
            args.warn_tol
        )
        placed.append(lig_i)
        print()

    # strip Rb and reorder core
    core = qd.copy()
    del core[[i for i,a in enumerate(core) if a.symbol=='Rb']]
    ordered_core = Atoms()
    for sym in ('Cs','Pb','Br','Cl'):
        idxs = [j for j,a in enumerate(core) if a.symbol==sym]
        if idxs:
            ordered_core += core[idxs]

    # merge & write
    final = ordered_core.copy()
    for lig in placed:
        final += lig
    io.write(args.out, final)
    print(f"\nWrote {args.out}  (Cs→Pb→Br→Cl → {len(placed)} ligands)")

if __name__=='__main__':
    main()

