
from ase.io import read, write
from xtb.ase.calculator import XTB
from sella import Sella
import traceback, os
inp = os.environ['TS_INPUT']
out = os.environ['TS_OUTPUT']
log = os.environ['TS_LOG']
traj = os.environ['TS_TRAJ']
atoms = read(inp)
atoms.calc = XTB(method='GFN2-xTB')
try:
    dyn = Sella(atoms, trajectory=traj, logfile=log, order=1)
    dyn.run(fmax=0.05, steps=250)
    write(out, atoms)
    print('FINAL_ENERGY', atoms.get_potential_energy())
except Exception as e:
    traceback.print_exc()
    write(out, atoms)
    raise
