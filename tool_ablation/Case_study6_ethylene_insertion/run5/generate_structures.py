import numpy as np
from math import pi, cos, sin, sqrt


def normalize(v):
    v = np.array(v, dtype=float)
    n = np.linalg.norm(v)
    if n < 1e-12:
        raise ValueError('zero vector')
    return v / n


def orthonormal_basis(n):
    n = normalize(n)
    tmp = np.array([0.0, 0.0, 1.0])
    if abs(np.dot(n, tmp)) > 0.9:
        tmp = np.array([0.0, 1.0, 0.0])
    u = normalize(np.cross(n, tmp))
    v = normalize(np.cross(n, u))
    return u, v


def cp_ring(centroid, normal, phase=0.0, R=1.19, CH=1.08):
    centroid = np.array(centroid, float)
    u, v = orthonormal_basis(normal)
    atoms = []
    for i in range(5):
        th = phase + 2*pi*i/5
        radial = cos(th)*u + sin(th)*v
        c = centroid + R*radial
        h = c + CH*radial
        atoms.append(("C", c))
        atoms.append(("H", h))
    return atoms


def add_ch3(carbon, neighbor_dirs, bond=1.09, phase=0.0):
    carbon = np.array(carbon, float)
    # neighbor_dirs: list of vectors from carbon to interacting heavy atoms
    axis = -normalize(np.sum([normalize(v) for v in neighbor_dirs], axis=0))
    e1, e2 = orthonormal_basis(axis)
    # tetrahedral angle with axis opposite heavy-atom direction group
    c = 1.0/3.0
    s = sqrt(1.0 - c*c)
    Hs = []
    for k in range(3):
        ph = phase + 2*pi*k/3
        d = normalize(c*axis + s*(cos(ph)*e1 + sin(ph)*e2))
        Hs.append(("H", carbon + bond*d))
    return Hs


def add_ch2_sp3(carbon, n1, n2, bond=1.09):
    carbon = np.array(carbon, float)
    u1 = normalize(n1)
    u2 = normalize(n2)
    cross = np.cross(u1, u2)
    if np.linalg.norm(cross) < 1e-6:
        # fallback
        axis = -u1
        e1, _ = orthonormal_basis(axis)
        H1 = normalize(axis + 0.95*e1)
        H2 = normalize(axis - 0.95*e1)
    else:
        p = normalize(cross)
        m = normalize(u1 + u2)
        c = (1.0/3.0) / max(np.dot(m, u1), 0.35)
        c = min(c, 0.95)
        s = sqrt(max(1.0 - c*c, 1e-8))
        H1 = normalize(-c*m + s*p)
        H2 = normalize(-c*m - s*p)
    return [("H", carbon + bond*H1), ("H", carbon + bond*H2)]


def add_ch2_sp2(carbon, carbon_neighbor, plane_normal, bond=1.09):
    carbon = np.array(carbon, float)
    u = normalize(carbon_neighbor)
    n = normalize(plane_normal)
    p = normalize(np.cross(n, u))
    H1 = normalize(-0.5*u + (sqrt(3)/2.0)*p)
    H2 = normalize(-0.5*u - (sqrt(3)/2.0)*p)
    return [("H", carbon + bond*H1), ("H", carbon + bond*H2)]


def write_xyz(path, atoms, comment=""):
    with open(path, 'w') as f:
        f.write(f"{len(atoms)}\n")
        f.write(comment + "\n")
        for el, xyz in atoms:
            x, y, z = xyz
            f.write(f"{el:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")


def dist(a, b):
    return float(np.linalg.norm(np.array(a) - np.array(b)))


def build_common_cp():
    zr = np.array([0.0, 0.0, 0.0])
    c1 = np.array([-1.45, 0.00, 1.75])
    c2 = np.array([-1.45, 0.00, -1.75])
    ring1 = cp_ring(c1, zr - c1, phase=0.35)
    ring2 = cp_ring(c2, zr - c2, phase=-0.25)
    return zr, ring1, ring2


def build_reactant():
    zr, ring1, ring2 = build_common_cp()
    atoms = [("Zr", zr)] + ring1 + ring2

    cm = np.array([1.30, -1.75, -0.45])
    ca = np.array([2.30,  0.70,  0.67])
    cb = np.array([2.30,  0.70, -0.67])

    atoms += [("C", cm), ("C", ca), ("C", cb)]
    atoms += add_ch3(cm, [zr - cm], phase=0.1)
    # Ethylene treated as sp2 in plane defined by C=C and Zr vector; use normal from cross(CC, Zr-centroid)
    cc = cb - ca
    plane_n = np.cross(cc, zr - (ca + cb)/2)
    if np.linalg.norm(plane_n) < 1e-6:
        plane_n = np.array([0.0, 1.0, 0.0])
    atoms += add_ch2_sp2(ca, cb - ca, plane_n)
    atoms += add_ch2_sp2(cb, ca - cb, plane_n)

    comment = (
        "Idealized reactant: [Cp2Zr(CH3)(eta2-C2H4)]+ ; Cp = unsubstituted cyclopentadienyl"
    )
    return atoms, comment, {
        'Zr-Cm': dist(zr, cm),
        'C=C': dist(ca, cb),
        'Zr-Ca': dist(zr, ca),
        'Zr-Cb': dist(zr, cb),
        'Cm...Cb': dist(cm, cb),
    }


def build_ts_like():
    zr, ring1, ring2 = build_common_cp()
    atoms = [("Zr", zr)] + ring1 + ring2

    ca = np.array([2.15,  0.45,  0.55])   # alpha carbon: remains metal-bound in product
    cb = np.array([3.08,  1.22, -0.04])   # beta carbon: receives migrating methyl
    cm = np.array([2.38, -0.52, -0.58])   # migrating methyl carbon

    atoms += [("C", ca), ("C", cb), ("C", cm)]
    # TS-like H placement: approximate only
    atoms += add_ch2_sp3(ca, zr - ca, cb - ca)
    atoms += add_ch2_sp3(cb, ca - cb, cm - cb)
    atoms += add_ch3(cm, [zr - cm, cb - cm], phase=0.2)

    comment = (
        "Idealized TS-like geometry for ethylene insertion into Zr-CH3; partial bonds implied by distances, not a computed TS"
    )
    return atoms, comment, {
        'Zr-Ca': dist(zr, ca),
        'Ca-Cb': dist(ca, cb),
        'Cb...Cm': dist(cb, cm),
        'Zr...Cm': dist(zr, cm),
    }


def build_product():
    zr, ring1, ring2 = build_common_cp()
    atoms = [("Zr", zr)] + ring1 + ring2

    ca = np.array([2.20, 0.40, 0.60])
    cb = np.array([3.05, 1.35, -0.25])
    cm = np.array([3.30, 2.45, -1.30])

    atoms += [("C", ca), ("C", cb), ("C", cm)]
    atoms += add_ch2_sp3(ca, zr - ca, cb - ca)
    atoms += add_ch2_sp3(cb, ca - cb, cm - cb)
    atoms += add_ch3(cm, [cb - cm], phase=0.3)

    comment = (
        "Idealized product: [Cp2Zr(CH2CH2CH3)]+ after 1,2-insertion of ethylene into the Zr-CH3 bond"
    )
    return atoms, comment, {
        'Zr-Ca': dist(zr, ca),
        'Ca-Cb': dist(ca, cb),
        'Cb-Cm': dist(cb, cm),
    }


def write_multi_xyz(path, entries):
    with open(path, 'w') as f:
        for atoms, comment, _ in entries:
            f.write(f"{len(atoms)}\n{comment}\n")
            for el, xyz in atoms:
                x, y, z = xyz
                f.write(f"{el:2s} {x:12.6f} {y:12.6f} {z:12.6f}\n")


def main():
    react = build_reactant()
    ts = build_ts_like()
    prod = build_product()
    write_xyz('reactant.xyz', react[0], react[1])
    write_xyz('ts_like.xyz', ts[0], ts[1])
    write_xyz('product.xyz', prod[0], prod[1])
    write_multi_xyz('zr_insertion_structures.xyz', [react, ts, prod])

    for name, data in [('reactant', react), ('ts_like', ts), ('product', prod)]:
        print(name)
        for k, v in data[2].items():
            print(f'  {k}: {v:.3f} Å')

if __name__ == '__main__':
    main()
