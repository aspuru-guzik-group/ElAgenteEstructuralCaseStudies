#!/usr/bin/env python3
"""TS optimization using Sella + xTB"""

import sys
import os
from ase.io import read, write
from ase.calculators.calculator import Calculator
import subprocess
import numpy as np

class XTBCalculator(Calculator):
    """Simple xTB calculator for ASE"""
    implemented_properties = ['energy', 'forces']
    
    def __init__(self, method='gfn2', **kwargs):
        Calculator.__init__(self, **kwargs)
        self.method = method
        
    def calculate(self, atoms=None, properties=['energy'], system_changes=None):
        Calculator.calculate(self, atoms, properties, system_changes)
        
        # Write atoms to temporary xyz file
        write('tmp_xtb.xyz', atoms)
        
        # Run xTB
        result = subprocess.run(
            ['xtb', 'tmp_xtb.xyz', '--' + self.method, '--grad'],
            capture_output=True, text=True
        )
        
        # Parse energy from output
        for line in result.stdout.split('\n'):
            if 'TOTAL ENERGY' in line:
                energy = float(line.split()[3]) * 27.211386245988  # Eh to eV
                
        # Parse forces from gradient file
        forces = []
        if os.path.exists('gradient'):
            with open('gradient', 'r') as f:
                lines = f.readlines()
            n_atoms = len(atoms)
            grad_start = None
            for i, line in enumerate(lines):
                if 'cycle' in line:
                    grad_start = i + n_atoms + 1
                    break
            if grad_start:
                for i in range(n_atoms):
                    parts = lines[grad_start + i].split()
                    grad = [float(x.replace('D', 'E')) for x in parts]
                    forces.append([-g * 27.211386245988 / 0.529177 for g in grad])  # Eh/Bohr to eV/Ang
        
        self.results['energy'] = energy
        self.results['forces'] = np.array(forces)

# Run TS optimization
from sella import Sella

input_file = sys.argv[1]
output_file = sys.argv[2]
log_file = sys.argv[3]

atoms = read(input_file)
atoms.calc = XTBCalculator(method='gfn2')

# Sella TS optimization
opt = Sella(atoms, logfile=log_file, trajectory=output_file.replace('.xyz', '.traj'))

try:
    converged = opt.run(fmax=0.02, steps=200)
    print(f"Converged: {converged}")
except Exception as e:
    print(f"Error during optimization: {e}")

# Save final structure
write(output_file, atoms)
print(f"Final energy: {atoms.get_potential_energy()} eV")
print(f"Structure saved to {output_file}")
