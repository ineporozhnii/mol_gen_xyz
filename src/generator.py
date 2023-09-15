import random
import numpy as np
import networkx as nx
from typing import Union
from sklearn.neighbors import KDTree
from src.utils import atomic_numbers_to_symbols

class Molecule:
    """
    Class to store properties of generated molecules
    """
    def __init__(self, 
                 atomic_numbers: list[int],
                 positions: np.ndarray
                 ):
        
        self.atomic_numbers = atomic_numbers
        self.positions = positions
        self.atomic_symbols = [atomic_numbers_to_symbols[number] for number in atomic_numbers]
        self.xyz_block = self._get_xyz_block()

    def _get_xyz_block(self):
        xyz_block = f"{len(self.atomic_numbers)}\n \n"
        for i, (symbol, position) in enumerate(zip(self.atomic_symbols, self.positions)):
            if i == len(self.atomic_numbers)-1:
                xyz_block += f"{symbol}    {position[0]:.5f}     {position[1]:.5f}     {position[2]:.5f}"
            else:
                xyz_block += f"{symbol}    {position[0]:.5f}     {position[1]:.5f}     {position[2]:.5f}\n"
        return xyz_block
    


def generate_mol(n_atoms_min: int,
                 n_atoms_max: int,
                 atomic_numbers: list[int],
                 unit_cell: list[float], 
                 min_distance: float,
                 max_distance: float,
                 no_disconnected_mols: True,
                 limit: int = 1000
                 ) -> Union[Molecule, bool]:
    """
    Finction to generate molecules with random atomic positions

    Args:
        n_atoms_min - minimum number of atoms in a molecule
        n_atoms_max - maximum number of atoms in a molecule
        atomic_numbers - allowed atomic numbers to choose from
        unit_cell - size of the unit cel
        min_distance - minimum allowed distance between atoms
        max_distance - maximum allowed distance between atoms
        no_disconnected_mols - if True, accepts only molecules with one subgraph 
        limit - maximum number of attempts to place atom in the unit cell

    Output:
        Molecule - generated molecule 
        accept - whether the generated molecule satisfies all conditions
    """

    accept = True
    if n_atoms_min == n_atoms_max:
        n_atoms = n_atoms_max
    else:
        n_atoms = np.random.randint(n_atoms_min, n_atoms_max+1)

    atoms = random.choices(atomic_numbers, k=n_atoms)

    check_coords = []
    system_coords = []

    coords = [unit_cell[0]/2, unit_cell[1]/2, unit_cell[2]/2]  #first atom in the center
    system_coords.append(coords)
    check_coords.append(coords)
    
    n_failed_atoms = 0
    while len(system_coords) < n_atoms:
        x = random.uniform(0, unit_cell[0])
        y = random.uniform(0, unit_cell[1])
        z = random.uniform(0, unit_cell[2])
        coords = [x, y, z]
        check_coords.append(coords)
        arr = np.array(check_coords)
        tree = KDTree(arr)

        if tree.query_radius(arr[-1].reshape(1, -1), r=min_distance, count_only=True) == 1:
            system_coords.append(coords)
        else:
            n_failed_atoms += 1
            del check_coords[-1]
        
        if n_failed_atoms > limit:
            return None, False

    system_coords = np.array(system_coords)

    bond_index = [[],[]]
    for idx in range(len(system_coords)):
        ind = []
        ind_radius = tree.query_radius(system_coords[idx].reshape(1, -1), r=max_distance)

        if len(ind_radius[0]) == 1:
            accept = False
        
        for j in range(len(ind_radius[0])):
            ind.append(ind_radius[0][j])

        for k in range(len(ind)):
            bond_index[0].append(idx)
            bond_index[1].append(ind[k])

    bonds = list(set([bond for bond in list(zip(bond_index[0], bond_index[1])) if bond[0]!=bond[1]]))

    networkx_graph = nx.Graph(bonds)
    subgraphs_list = list(nx.connected_components(networkx_graph))

    if no_disconnected_mols and len(subgraphs_list) > 1:
        accept = False

    mol = Molecule(atomic_numbers = atoms, positions=system_coords)

    return mol, accept





    
