import sys, os, json, traceback
from ase.io import read, write
from xtb.ase.calculator import XTB
from sella import Sella

inp, out_xyz, log_path, traj_path = sys.argv[1:5]
atoms = read(inp)
atoms.calc = XTB(method='GFN2-xTB', charge=0, uhf=0)
try:
    opt = Sella(atoms, trajectory=traj_path, logfile=log_path, order=1, internal=True, allow_fragments=True)
    opt.run(fmax=0.05, steps=150)
    write(out_xyz, atoms)
    print('DONE')
except Exception as e:
    traceback.print_exc()
    write(out_xyz, atoms)
    sys.exit(1)
