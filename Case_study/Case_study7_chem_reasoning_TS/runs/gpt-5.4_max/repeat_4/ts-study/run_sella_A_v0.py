import json, traceback
from pathlib import Path
import numpy as np
from ase.io import read, write
from xtb.ase.calculator import XTB
from sella import Sella
root = Path(__file__).resolve().parent
base = root / 'calculations' / 'pathway_A'
start = base / 'ts_guess_constrained.xyz'
target = base / 'pseudo_product_like.xyz'
out = root / 'structures' / 'pathway_A_TS_optimized.xyz'
traj = base / 'sellaA_v0.traj'
log = base / 'sellaA_v0.log'
resj = base / 'sellaA_v0_result.json'
try:
    atoms = read(str(start))
    targ = read(str(target))
    v0 = (targ.get_positions() - atoms.get_positions()).reshape(-1)
    v0 = v0 / np.linalg.norm(v0)
    atoms.calc = XTB(method='GFN2-xTB')
    opt = Sella(
        atoms,
        trajectory=str(traj),
        logfile=str(log),
        order=1,
        internal=False,
        nsteps_per_diag=3,
        diag_every_n=1,
        eta=1e-4,
        allow_fragments=True,
        v0=v0,
    )
    opt.run(fmax=0.03, steps=150)
    write(str(out), atoms)
    e = atoms.get_potential_energy()
    with open(resj, 'w') as fh:
        json.dump({'status':'ok','energy_eV':e}, fh)
except Exception as e:
    with open(resj, 'w') as fh:
        json.dump({'status':'error','error':str(e), 'traceback': traceback.format_exc()}, fh)
    raise
