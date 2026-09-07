
import numpy as np
from math import cos, sin, pi
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt


def norm(v):
    v = np.array(v, dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-10:
        raise ValueError('zero vector')
    return v/n


def orthonormal_basis_from_normal(n):
    n = norm(n)
    ref = np.array([0.0, 1.0, 0.0])
    if abs(np.dot(n, ref)) > 0.9:
        ref = np.array([1.0, 0.0, 0.0])
    e1 = norm(np.cross(n, ref))
    e2 = norm(np.cross(n, e1))
    return e1, e2


def place_ch3_hs(cpos, attached_to, ch=1.09, phase=0.0):
    cpos = np.array(cpos, float)
    u = norm(np.array(attached_to, float) - cpos)
    e1, e2 = orthonormal_basis_from_normal(u)
    hs = []
    for ang in [phase, phase + 2*pi/3, phase + 4*pi/3]:
        d = -u/3.0 + (2*np.sqrt(2)/3.0)*(cos(ang)*e1 + sin(ang)*e2)
        hs.append(cpos + ch*norm(d))
    return hs


def best_ch3_hs(cpos, attached_to, avoid_points, ch=1.09, phases=72):
    best = None
    best_score = -1.0
    for phase in np.linspace(0, 2*pi, phases, endpoint=False):
        hs = place_ch3_hs(cpos, attached_to, ch=ch, phase=phase)
        score = min(np.linalg.norm(h - p) for h in hs for p in avoid_points)
        if score > best_score:
            best_score = score
            best = hs
    return best


def place_ch2_hs(cpos, neighbor1, neighbor2, ch=1.09, spread=0.95):
    cpos = np.array(cpos, float)
    u1 = norm(np.array(neighbor1, float) - cpos)
    u2 = norm(np.array(neighbor2, float) - cpos)
    b = -(u1 + u2)
    if np.linalg.norm(b) < 1e-6:
        ref = np.array([0.0, 0.0, 1.0])
        if abs(np.dot(u1, ref)) > 0.9:
            ref = np.array([0.0, 1.0, 0.0])
        b = np.cross(u1, ref)
    b = norm(b)
    n = np.cross(u1, u2)
    if np.linalg.norm(n) < 1e-6:
        e1, _ = orthonormal_basis_from_normal(b)
        n = e1
    n = norm(n)
    d1 = norm(b + spread*n)
    d2 = norm(b - spread*n)
    return [cpos + ch*d1, cpos + ch*d2]


def add_cp_ring(atoms, bonds, centroid, normal, phase=0.0, ring_radius=1.68, ch=1.08):
    e1, e2 = orthonormal_basis_from_normal(normal)
    carbons = []
    for i in range(5):
        ang = phase + 2*pi*i/5
        radial = cos(ang)*e1 + sin(ang)*e2
        c = centroid + ring_radius*radial
        atoms.append(('C', c))
        ci = len(atoms)-1
        carbons.append(ci)
        h = c + ch*radial
        atoms.append(('H', h))
        bonds.append((ci, len(atoms)-1, 'solid'))
    for i in range(5):
        bonds.append((carbons[i], carbons[(i+1)%5], 'solid'))


def build_common_scaffold():
    atoms=[('Zr', np.array([0.0,0.0,0.0]))]
    bonds=[]
    d_cp = 2.30
    ang = np.deg2rad(65.0)
    c1 = np.array([-d_cp*cos(ang), 0.0, d_cp*sin(ang)])
    c2 = np.array([-d_cp*cos(ang), 0.0,-d_cp*sin(ang)])
    add_cp_ring(atoms,bonds,c1,c1,phase=0.2)
    add_cp_ring(atoms,bonds,c2,c2,phase=-0.3)
    return atoms,bonds


def build_reactant():
    atoms,bonds = build_common_scaffold()
    zr = atoms[0][1]
    # methyl
    c_me = np.array([2.15,-0.65,-0.35])
    me_idx = len(atoms); atoms.append(('C', c_me)); bonds.append((0, me_idx, 'solid'))
    # delay H placement until alkene carbons are defined for clash-avoidance
    # coordinated ethylene: proximal carbon near Zr, alkene plane chosen to avoid Me clash
    c_a = np.array([2.35,0.70,0.45])
    c_b = np.array([3.63,1.05,0.45])
    a_idx = len(atoms); atoms.append(('C', c_a))
    b_idx = len(atoms); atoms.append(('C', c_b))
    bonds.append((a_idx, b_idx, 'double'))
    bonds.append((0, a_idx, 'partial'))
    u_ab = norm(c_b - c_a)
    plane_n = np.array([0.0,1.0,0.0])
    v = norm(np.cross(plane_n, u_ab))
    for cpos, cidx, u in [(c_a, a_idx, u_ab), (c_b, b_idx, -u_ab)]:
        for sgn in [+1,-1]:
            d = norm(-0.5*u + sgn*(np.sqrt(3)/2.0)*v)
            atoms.append(('H', cpos + 1.08*d)); bonds.append((cidx, len(atoms)-1, 'solid'))
    avoid = [p for _,p in atoms if True] + [c_a, c_b]
    for h in best_ch3_hs(c_me, zr, avoid):
        atoms.append(('H', h)); bonds.append((me_idx, len(atoms)-1, 'solid'))
    return atoms,bonds


def build_ts_like():
    atoms,bonds = build_common_scaffold()
    zr = atoms[0][1]
    c_a = np.array([2.20,0.55,0.45])
    c_b = np.array([3.55,1.10,0.20])
    c_me = np.array([2.35,-0.35,-0.45])
    a_idx = len(atoms); atoms.append(('C', c_a))
    b_idx = len(atoms); atoms.append(('C', c_b))
    me_idx = len(atoms); atoms.append(('C', c_me))
    bonds += [(0,a_idx,'partial'), (0,me_idx,'partial'), (a_idx,b_idx,'partial2'), (b_idx,me_idx,'partial')]
    for h in place_ch2_hs(c_a, zr, c_b):
        atoms.append(('H', h)); bonds.append((a_idx, len(atoms)-1, 'solid'))
    for h in place_ch2_hs(c_b, c_a, c_me):
        atoms.append(('H', h)); bonds.append((b_idx, len(atoms)-1, 'solid'))
    avoid = [p for _,p in atoms if True] + [c_a, c_b]
    # orient CH3 H atoms away from the insertion center
    for h in best_ch3_hs(c_me, (zr + c_b)/2.0, avoid):
        atoms.append(('H', h)); bonds.append((me_idx, len(atoms)-1, 'solid'))
    return atoms,bonds


def build_product():
    atoms,bonds = build_common_scaffold()
    zr = atoms[0][1]
    c1 = np.array([2.15,0.45,0.45])
    c2 = np.array([3.55,1.05,0.15])
    c3 = np.array([4.10,2.30,-0.55])
    i1 = len(atoms); atoms.append(('C', c1)); bonds.append((0,i1,'solid'))
    i2 = len(atoms); atoms.append(('C', c2)); bonds.append((i1,i2,'solid'))
    i3 = len(atoms); atoms.append(('C', c3)); bonds.append((i2,i3,'solid'))
    for h in place_ch2_hs(c1, zr, c2):
        atoms.append(('H', h)); bonds.append((i1, len(atoms)-1, 'solid'))
    for h in place_ch2_hs(c2, c1, c3):
        atoms.append(('H', h)); bonds.append((i2, len(atoms)-1, 'solid'))
    avoid = [p for _,p in atoms if True]
    for h in best_ch3_hs(c3, c2, avoid):
        atoms.append(('H', h)); bonds.append((i3, len(atoms)-1, 'solid'))
    return atoms,bonds


def write_xyz(fname, atoms, comment):
    with open(fname,'w') as f:
        f.write(f"{len(atoms)}\n{comment}\n")
        for el,p in atoms:
            f.write(f"{el:2s} {p[0]:12.6f} {p[1]:12.6f} {p[2]:12.6f}\n")


def pairwise_min_dist(atoms):
    coords=np.array([p for _,p in atoms])
    mind=1e9; pair=None
    for i in range(len(coords)):
        for j in range(i+1,len(coords)):
            d=np.linalg.norm(coords[i]-coords[j])
            if d<mind:
                mind=d; pair=(i,j)
    return mind,pair


def draw_structure(ax, atoms, bonds, title):
    colors={'C':'#444444','H':'#d9d9d9','Zr':'#5b8ff9'}
    sizes={'C':80,'H':30,'Zr':220}
    coords=np.array([p for _,p in atoms])
    for a,b,style in bonds:
        pa=atoms[a][1]; pb=atoms[b][1]
        if style=='solid':
            ax.plot([pa[0],pb[0]],[pa[1],pb[1]],[pa[2],pb[2]],color='k',lw=1.6)
        elif style=='double':
            ax.plot([pa[0],pb[0]],[pa[1],pb[1]],[pa[2],pb[2]],color='k',lw=2.4)
        elif style=='partial':
            ax.plot([pa[0],pb[0]],[pa[1],pb[1]],[pa[2],pb[2]],color='k',lw=1.4,ls='--')
        elif style=='partial2':
            ax.plot([pa[0],pb[0]],[pa[1],pb[1]],[pa[2],pb[2]],color='dimgray',lw=1.4,ls='-.')
    for el in ['Zr','C','H']:
        idx=[i for i,(e,_) in enumerate(atoms) if e==el]
        if idx:
            xyz=coords[idx]
            ax.scatter(xyz[:,0],xyz[:,1],xyz[:,2],s=sizes[el],c=colors[el],edgecolors='k',linewidths=0.4,depthshade=True)
    ax.set_title(title,fontsize=10)
    ax.view_init(elev=18, azim=-62)
    maxr = np.max(np.ptp(coords,axis=0))/2 + 0.8
    ctr = coords.mean(axis=0)
    ax.set_xlim(ctr[0]-maxr, ctr[0]+maxr)
    ax.set_ylim(ctr[1]-maxr, ctr[1]+maxr)
    ax.set_zlim(ctr[2]-maxr, ctr[2]+maxr)
    ax.set_axis_off()


def main():
    data = {
        'reactant.xyz': (build_reactant(), '[Cp2Zr(CH3)(ethylene)]+ approximate 3D model; counterion omitted'),
        'ts_like.xyz': (build_ts_like(), 'Approximate migratory-insertion TS-like geometry for ethylene insertion into Zr-CH3; not an optimized TS'),
        'product.xyz': (build_product(), '[Cp2Zr(n-propyl)]+ approximate 3D model; counterion omitted'),
    }
    fig=plt.figure(figsize=(12,4))
    lines=[]
    for i,(fname,((atoms,bonds),comment)) in enumerate(data.items(), start=1):
        write_xyz(fname, atoms, comment)
        md,pair = pairwise_min_dist(atoms)
        lines.append(f'{fname}: atoms={len(atoms)}, min interatomic distance={md:.2f} Å between indices {pair}')
        ax=fig.add_subplot(1,3,i,projection='3d')
        draw_structure(ax, atoms, bonds, fname[:-4])
    plt.tight_layout()
    plt.savefig('generated_structures.png', dpi=200, bbox_inches='tight')
    with open('structure_generation_report.txt','w') as f:
        f.write('\n'.join(lines)+'\n')
        f.write('These are approximate hand-built geometries intended as plausible 3D starting structures, not optimized stationary points.\n')

if __name__ == '__main__':
    main()
