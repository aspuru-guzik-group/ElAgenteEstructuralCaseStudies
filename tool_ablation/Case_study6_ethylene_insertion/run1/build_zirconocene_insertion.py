import math
import numpy as np
from pathlib import Path

# ---------- vector helpers ----------
def norm(v):
    v = np.array(v, dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-12:
        raise ValueError('zero vector')
    return v / n

def orthonormal_basis(n):
    n = norm(n)
    ref = np.array([1.0, 0.0, 0.0])
    if abs(np.dot(n, ref)) > 0.85:
        ref = np.array([0.0, 1.0, 0.0])
    u = norm(np.cross(n, ref))
    v = norm(np.cross(n, u))
    return u, v

# ---------- geometry builders ----------
def cp_ring(centroid, metal=np.zeros(3), c_r=1.19, ch=1.08, phase_deg=0.0, ring_name='Cp'):
    centroid = np.array(centroid, float)
    n = norm(metal - centroid)
    u, v = orthonormal_basis(n)
    phase = math.radians(phase_deg)
    atoms = []
    carbons = []
    for i in range(5):
        ang = phase + 2 * math.pi * i / 5
        radial = math.cos(ang) * u + math.sin(ang) * v
        c = centroid + c_r * radial
        carbons.append(c)
        atoms.append(('C', c, f'{ring_name}C{i+1}'))
    for i, c in enumerate(carbons):
        radial = norm(c - centroid)
        h = centroid + (c_r + ch) * radial
        atoms.append(('H', h, f'{ring_name}H{i+1}'))
    return atoms

def add_sp3_one_bond(center, attached, bond_len=1.09, phase_deg=0.0, label_prefix='H'):
    center = np.array(center, float)
    n = norm(np.array(attached, float) - center)
    u, v = orthonormal_basis(n)
    atoms = []
    phase = math.radians(phase_deg)
    for i, phi_deg in enumerate((0, 120, 240), start=1):
        phi = phase + math.radians(phi_deg)
        d = (-1.0/3.0) * n + (2.0 * math.sqrt(2.0) / 3.0) * (math.cos(phi) * u + math.sin(phi) * v)
        atoms.append(('H', center + bond_len * norm(d), f'{label_prefix}{i}'))
    return atoms

def add_sp3_two_bonds(center, attached1, attached2, bond_len=1.09, twist_deg=0.0, label_prefix='H'):
    center = np.array(center, float)
    n1 = norm(np.array(attached1, float) - center)
    n2 = norm(np.array(attached2, float) - center)
    s = n1 + n2
    if np.linalg.norm(s) < 1e-8:
        b = np.array([1.0, 0.0, 0.0])
    else:
        b = norm(s)
    p = np.cross(n1, n2)
    if np.linalg.norm(p) < 1e-8:
        p = orthonormal_basis(b)[0]
    else:
        p = norm(p)
    q = norm(np.cross(p, b))
    # base tetrahedral-like directions around the bisector
    # each H should be mostly opposite the bisector with a perpendicular component
    a = 0.58  # ~1/sqrt(3), suitable for tetrahedral-like CH2 placement
    c = math.sqrt(max(1e-8, 1.0 - a*a))
    tw = math.radians(twist_deg)
    perp1 = math.cos(tw) * p + math.sin(tw) * q
    perp2 = -math.cos(tw) * p + math.sin(tw) * q
    d1 = norm(-a * b + c * perp1)
    d2 = norm(-a * b + c * perp2)
    return [('H', center + bond_len * d1, f'{label_prefix}1'), ('H', center + bond_len * d2, f'{label_prefix}2')]

def add_sp2_alkene_hydrogens(center, other_c, plane_normal, bond_len=1.08, label_prefix='H'):
    center = np.array(center, float)
    other_c = np.array(other_c, float)
    t = norm(other_c - center)
    n = norm(plane_normal)
    w = np.cross(n, t)
    if np.linalg.norm(w) < 1e-8:
        w = orthonormal_basis(t)[0]
    else:
        w = norm(w)
    d1 = norm(-0.5 * t + (math.sqrt(3)/2.0) * w)
    d2 = norm(-0.5 * t - (math.sqrt(3)/2.0) * w)
    return [('H', center + bond_len * d1, f'{label_prefix}1'), ('H', center + bond_len * d2, f'{label_prefix}2')]

# ---------- writing and checks ----------
def write_xyz(filename, atoms, title):
    with open(filename, 'w') as f:
        f.write(f"{len(atoms)}\n{title}\n")
        for el, xyz, label in atoms:
            x, y, z = xyz
            f.write(f"{el:2s} {x:12.6f} {y:12.6f} {z:12.6f}  # {label}\n")

def bondset_cp(ring_prefix):
    pairs = set()
    for i in range(1,6):
        j = i+1 if i < 5 else 1
        pairs.add(tuple(sorted((f'{ring_prefix}C{i}', f'{ring_prefix}C{j}'))))
        pairs.add(tuple(sorted((f'{ring_prefix}C{i}', f'{ring_prefix}H{i}'))))
    return pairs

def find_nonbonded_clash(atoms, bonded_pairs):
    coords = {label: np.array(xyz) for el, xyz, label in atoms}
    elems = {label: el for el, xyz, label in atoms}
    labels = list(coords)
    worst = None
    for i in range(len(labels)):
        for j in range(i+1, len(labels)):
            li, lj = labels[i], labels[j]
            if tuple(sorted((li, lj))) in bonded_pairs:
                continue
            d = np.linalg.norm(coords[li] - coords[lj])
            cutoff = 1.05 if 'H' in (elems[li], elems[lj]) else 1.25
            if d < cutoff:
                if worst is None or d < worst[0]:
                    worst = (d, li, lj)
    return worst

# ---------- base framework ----------
Zr = np.array([0.0, 0.0, 0.0])
cp1_cent = np.array([-1.20,  1.80,  0.80])
cp2_cent = np.array([-1.20, -1.80, -0.80])
cp_atoms = cp_ring(cp1_cent, metal=Zr, phase_deg=18, ring_name='CpA') + cp_ring(cp2_cent, metal=Zr, phase_deg=-8, ring_name='CpB')
cp_bonds = bondset_cp('CpA') | bondset_cp('CpB')

# ---------- reactant ----------
reactant = [('Zr', Zr, 'Zr')] + cp_atoms
Cm = np.array([2.05, -0.75, 0.85])
Ca = np.array([2.00, 1.10, -0.35])
Cb = np.array([2.55, 0.20, -1.18])
reactant += [('C', Cm, 'MeC'), ('C', Ca, 'AlkeneC1'), ('C', Cb, 'AlkeneC2')]
reactant += add_sp3_one_bond(Cm, Zr, phase_deg=25, label_prefix='MeH')
plane_n = np.cross(Cb - Ca, Zr - (Ca + Cb) / 2.0)
reactant += add_sp2_alkene_hydrogens(Ca, Cb, plane_n, label_prefix='C1H')
reactant += add_sp2_alkene_hydrogens(Cb, Ca, plane_n, label_prefix='C2H')
react_bonds = set(cp_bonds)
react_bonds |= {tuple(sorted(x)) for x in [('Zr','MeC'), ('AlkeneC1','AlkeneC2'), ('MeC','MeH1'), ('MeC','MeH2'), ('MeC','MeH3'), ('AlkeneC1','C1H1'), ('AlkeneC1','C1H2'), ('AlkeneC2','C2H1'), ('AlkeneC2','C2H2')]}

# ---------- TS-like ----------
TS = [('Zr', Zr, 'Zr')] + cp_atoms
Cm_ts = np.array([2.20, -0.70, 0.30])
Ca_ts = np.array([2.15, 0.85, -0.45])
Cb_ts = np.array([3.15, 0.05, -1.15])
TS += [('C', Cm_ts, 'MigratingMeC'), ('C', Ca_ts, 'AlphaC'), ('C', Cb_ts, 'BetaC')]
TS += add_sp3_one_bond(Cm_ts, Zr, phase_deg=0, label_prefix='MeH')
TS += add_sp3_two_bonds(Ca_ts, Zr, Cb_ts, twist_deg=-50, label_prefix='AlphaH')
TS += add_sp3_two_bonds(Cb_ts, Ca_ts, Cm_ts, twist_deg=-50, label_prefix='BetaH')
TS_bonds = set(cp_bonds)
TS_bonds |= {tuple(sorted(x)) for x in [
    ('Zr','MigratingMeC'), ('Zr','AlphaC'), ('AlphaC','BetaC'), ('MigratingMeC','BetaC'),
    ('MigratingMeC','MeH1'), ('MigratingMeC','MeH2'), ('MigratingMeC','MeH3'),
    ('AlphaC','AlphaH1'), ('AlphaC','AlphaH2'), ('BetaC','BetaH1'), ('BetaC','BetaH2')
]}

# ---------- product ----------
product = [('Zr', Zr, 'Zr')] + cp_atoms
C1 = np.array([1.95, 1.05, -0.45])
d12 = norm(np.array([0.74, -0.62, 0.17]))
C2 = C1 + 1.53 * d12
# choose a bent C2-C3 direction giving a realistic chain angle
# approximately 110 deg from C2->C1 and pointing away from the metal pocket
v23 = norm(np.array([0.70, -0.10, -0.70]))
C3 = C2 + 1.53 * v23
product += [('C', C1, 'PropylC1'), ('C', C2, 'PropylC2'), ('C', C3, 'PropylC3')]
product += add_sp3_two_bonds(C1, Zr, C2, twist_deg=20, label_prefix='C1H')
product += add_sp3_two_bonds(C2, C1, C3, twist_deg=-25, label_prefix='C2H')
product += add_sp3_one_bond(C3, C2, phase_deg=40, label_prefix='C3H')
prod_bonds = set(cp_bonds)
prod_bonds |= {tuple(sorted(x)) for x in [
    ('Zr','PropylC1'), ('PropylC1','PropylC2'), ('PropylC2','PropylC3'),
    ('PropylC1','C1H1'), ('PropylC1','C1H2'), ('PropylC2','C2H1'), ('PropylC2','C2H2'),
    ('PropylC3','C3H1'), ('PropylC3','C3H2'), ('PropylC3','C3H3')
]}

# ---------- write ----------
write_xyz('reactant.xyz', reactant, '[Cp2Zr(CH3)(eta2-C2H4)]+ approximate explicit-atom 3D model')
write_xyz('ts_like.xyz', TS, 'TS-like [Cp2Zr(CH3)(C2H4)]+ insertion geometry; approximate, not a validated TS')
write_xyz('product.xyz', product, '[Cp2Zr(CH2CH2CH3)]+ approximate explicit-atom 3D model')

# ---------- checks ----------
def distance(a,b):
    return float(np.linalg.norm(a-b))

def coords_map(atoms):
    return {label: np.array(xyz) for el, xyz, label in atoms}

for name, atoms, bonds in [('reactant', reactant, react_bonds), ('ts_like', TS, TS_bonds), ('product', product, prod_bonds)]:
    m = coords_map(atoms)
    print(f'[{name}]')
    if name == 'reactant':
        print('  Zr-MeC      ', round(distance(m['Zr'], m['MeC']), 3))
        print('  Zr-C1       ', round(distance(m['Zr'], m['AlkeneC1']), 3))
        print('  Zr-C2       ', round(distance(m['Zr'], m['AlkeneC2']), 3))
        print('  C1=C2       ', round(distance(m['AlkeneC1'], m['AlkeneC2']), 3))
    elif name == 'ts_like':
        print('  Zr-MeC      ', round(distance(m['Zr'], m['MigratingMeC']), 3))
        print('  Zr-Alpha    ', round(distance(m['Zr'], m['AlphaC']), 3))
        print('  Alpha-Beta  ', round(distance(m['AlphaC'], m['BetaC']), 3))
        print('  Me...Beta   ', round(distance(m['MigratingMeC'], m['BetaC']), 3))
    else:
        print('  Zr-C1       ', round(distance(m['Zr'], m['PropylC1']), 3))
        print('  C1-C2       ', round(distance(m['PropylC1'], m['PropylC2']), 3))
        print('  C2-C3       ', round(distance(m['PropylC2'], m['PropylC3']), 3))
    clash = find_nonbonded_clash(atoms, bonds)
    print('  clash       ', clash)

# ---------- previews ----------
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
colors = {'Zr':'#2b6cb0','C':'#1f2937','H':'#cbd5e1'}
sizes = {'Zr':220,'C':65,'H':24}

def draw(name, atoms, bonds):
    fig = plt.figure(figsize=(5,5))
    ax = fig.add_subplot(111, projection='3d')
    for el in ['Zr','C','H']:
        pts = np.array([xyz for e, xyz, _ in atoms if e == el])
        if len(pts):
            ax.scatter(pts[:,0], pts[:,1], pts[:,2], s=sizes[el], c=colors[el], depthshade=True)
    dct = {label: np.array(xyz) for _, xyz, label in atoms}
    for a,b in bonds:
        if a == 'Zr' and b.startswith('Cp'):
            continue
        xa, xb = dct[a], dct[b]
        ax.plot([xa[0],xb[0]],[xa[1],xb[1]],[xa[2],xb[2]], color='k', lw=1.3)
    # stylized Zr-Cp centroid contacts
    for cent in [cp1_cent, cp2_cent]:
        ax.plot([0, cent[0]], [0, cent[1]], [0, cent[2]], color='gray', lw=1.1, alpha=0.7)
    allpts = np.array([xyz for _, xyz, _ in atoms])
    mins = allpts.min(axis=0); maxs = allpts.max(axis=0)
    center = (mins + maxs) / 2.0
    span = max(maxs - mins) * 0.62
    ax.set_xlim(center[0]-span, center[0]+span)
    ax.set_ylim(center[1]-span, center[1]+span)
    ax.set_zlim(center[2]-span, center[2]+span)
    ax.set_axis_off()
    ax.view_init(elev=18, azim=-58)
    plt.tight_layout()
    plt.savefig(f'{name}.png', dpi=180, bbox_inches='tight', pad_inches=0.02)
    plt.close(fig)

for name, atoms, bonds in [('reactant', reactant, react_bonds), ('ts_like', TS, TS_bonds), ('product', product, prod_bonds)]:
    draw(name, atoms, bonds)

with open('STRUCTURE_NOTES.txt', 'w') as f:
    f.write('Files:\n')
    f.write('  reactant.xyz : [Cp2Zr(CH3)(eta2-C2H4)]+ approximate 3D structure\n')
    f.write('  ts_like.xyz  : insertion TS-like geometry (illustrative, not optimized)\n')
    f.write('  product.xyz  : [Cp2Zr(CH2CH2CH3)]+ approximate 3D structure\n\n')
    f.write('Interpretation: migratory insertion of coordinated ethylene into the zirconium-methyl bond of a cationic zirconocene complex.\n')
    f.write('The TS-like model is schematic/approximate and should be refined before any energetic use.\n')
