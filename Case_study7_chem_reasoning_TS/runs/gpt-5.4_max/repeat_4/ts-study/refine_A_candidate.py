import json, traceback
from pathlib import Path
from ase.io import read, write
from xtb.ase.calculator import XTB
from sella import Sella
root = Path(__file__).resolve().parent
base = root / 'calculations' / 'pathway_A'
start = base / 'test_ts.xyz'
out = root / 'structures' / 'pathway_A_TS_optimized_refined.xyz'
traj = base / 'refineA.traj'
log = base / 'refineA.log'
resj = base / 'refineA_result.json'
try:
    atoms = read(str(start))
    atoms.calc = XTB(method='GFN2-xTB')
    opt = Sella(atoms, trajectory=str(traj), logfile=str(log), order=1, internal=False, allow_fragments=True, diag_every_n=1)
    opt.run(fmax=0.02, steps=80)
    write(str(out), atoms)
    with open(resj,'w') as fh:
        json.dump({'status':'ok','energy_eV':atoms.get_potential_energy()},fh)
except Exception as e:
    with open(resj,'w') as fh:
        json.dump({'status':'error','error':str(e), 'traceback': traceback.format_exc()},fh)
    raise
