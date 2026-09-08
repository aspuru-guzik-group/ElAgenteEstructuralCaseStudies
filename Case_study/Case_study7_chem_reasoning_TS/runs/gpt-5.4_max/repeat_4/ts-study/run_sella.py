import sys, json, traceback
from pathlib import Path
from ase.io import read, write
from xtb.ase.calculator import XTB
from sella import Sella

label = sys.argv[1]
root = Path(__file__).resolve().parent
base = root / 'calculations' / f'pathway_{label}'
start = base / 'ts_guess_constrained.xyz'
out = root / 'structures' / f'pathway_{label}_TS_optimized.xyz'
traj = base / 'sella_opt.traj'
log = base / 'sella_opt.log'
resj = base / 'sella_result.json'
try:
    atoms = read(str(start))
    atoms.calc = XTB(method='GFN2-xTB')
    opt = Sella(
        atoms,
        trajectory=str(traj),
        logfile=str(log),
        order=1,
        internal=True,
        nsteps_per_diag=3,
        diag_every_n=2,
        eta=1e-4,
        allow_fragments=True,
    )
    opt.run(fmax=0.03, steps=200)
    write(str(out), atoms)
    e = atoms.get_potential_energy()
    with open(resj, 'w') as fh:
        json.dump({'status':'ok','energy_eV':e}, fh)
except Exception as e:
    with open(resj, 'w') as fh:
        json.dump({'status':'error','error':str(e), 'traceback': traceback.format_exc()}, fh)
    raise
