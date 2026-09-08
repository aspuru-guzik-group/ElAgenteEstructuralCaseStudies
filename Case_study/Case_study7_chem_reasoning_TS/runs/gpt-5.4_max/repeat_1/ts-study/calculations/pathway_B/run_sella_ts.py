from ase.io import read, write
from xtb.ase.calculator import XTB
from sella import Sella

atoms = read('../../structures/pathway_B_TS_Initial_Guess.xyz')
atoms.calc = XTB(method='GFN2-xTB')

dyn = Sella(atoms,
            logfile='sella.log',
            trajectory='sella.traj',
            order=1,
            internal=False,
            allow_fragments=True,
            eta=1e-4)
dyn.run(fmax=0.05, steps=150)
write('../../structures/pathway_B_TS_optimized.xyz', atoms)
print('Final energy eV', atoms.get_potential_energy())
