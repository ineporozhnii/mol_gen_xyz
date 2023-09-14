import random
import numpy as np

class Molecule:
    def __init__(self, 
                 atomic_numbers: list[int],
                 positions: np.ndarray
                 ):
        
        self.atomic_numbers = atomic_numbers
        self.positions = positions
    


def generate_mol(n_atoms_min: int,
                 n_atoms_max: int,
                 atomic_numbers: list[int]
                 ) -> Molecule:
    
    if n_atoms_min == n_atoms_max:
        n_atoms = n_atoms_max
    else:
        n_atoms = np.random.randint(n_atoms_min, n_atoms_max+1)

    atoms = random.choices(atomic_numbers, k=n_atoms)