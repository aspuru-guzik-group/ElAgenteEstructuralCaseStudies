import sys, json, traceback
from ase.io import read
from xtb.ase.calculator import XTB
from ase.vibrations import Vibrations

inp, out_json, workdir = sys.argv[1:4]
atoms = read(inp)
atoms.calc = XTB(method='GFN2-xTB', charge=0, uhf=0)
try:
    vib = Vibrations(atoms, name=workdir+'/vib')
    vib.run()
    freqs = vib.get_frequencies()
    energies = vib.get_energies()
    with open(out_json,'w') as f:
        json.dump({'frequencies_cm-1':[float(x) for x in freqs], 'energies_eV':[float(x.real) if hasattr(x,'real') else float(x) for x in energies]}, f, indent=2)
    vib.summary(log=workdir+'/vib_summary.txt')
    print('DONE')
except Exception:
    traceback.print_exc()
    sys.exit(1)
