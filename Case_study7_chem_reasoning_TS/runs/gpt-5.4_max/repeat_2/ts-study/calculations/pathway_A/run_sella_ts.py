from ase.io import read, write
from xtb.ase.calculator import XTB
from sella import Sella

ats = read('../../structures/pathway_A_TS_Initial_Guess.xyz')
ats.calc = XTB(method='GFN2-xTB')
opt = Sella(
    ats,
    logfile='sella.log',
    trajectory='sella.traj',
    order=1,
    internal=True,
    eta=1e-4,
    nsteps_per_diag=3,
    diag_every_n=1,
)
opt.run(fmax=0.05, steps=200)
write('../../structures/pathway_A_TS_optimized.xyz', ats)
print('Final energy eV', ats.get_potential_energy())
