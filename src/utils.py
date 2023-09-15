import warnings
import py3Dmol
import numpy as np
from rdkit import Chem
from typing import Union
from stmol import showmol
from ase import Atoms
from ase.optimize import BFGS
from ase.calculators.lj import LennardJones
from rdkit.Chem import rdDetermineBonds
from rdkit.Chem.Draw import IPythonConsole
IPythonConsole.ipython_3d = True


def check_input(config_dict: dict) -> dict[str, Union[int, float, list, str]]:
    clean_config_dict = {}
    allowed_arguments = get_allowed_arguments()
    
    for key, value in config_dict.items():
        if key in allowed_arguments.keys():
            if not isinstance(value, allowed_arguments[key]):
                 warnings.warn(f"Invalid config argument for {key} will be ignored. Expected: {allowed_arguments[key]}, got: {type(value)}", stacklevel=3)
            else:
                clean_config_dict[key] = value
        else:
            warnings.warn(f"Unrecognized config argument {key} will be ignored.", stacklevel=3)
    
    return clean_config_dict
        

def get_allowed_arguments() -> dict[str, type]:
    return {"n_mols": int, "n_atoms_min": int, "n_atoms_max": int, "atomic_numbers": list, "unit_cell": list, "min_distance": float, "max_distance": float, "random_seed": int, "name_id": str, "no_disconnected_mols": bool}


def get_connected_mol(xyz_block: str):
    print(xyz_block)
    mol = Chem.rdmolfiles.MolFromXYZBlock(xyz_block)
    connected_mol = Chem.Mol(mol)
    rdDetermineBonds.DetermineConnectivity(connected_mol)
    print(type(connected_mol))
    return connected_mol


def display_molecule(mol_block):
    xyzview = py3Dmol.view(width=400,height=400)
    IPythonConsole.addMolToView(mol_block, xyzview)
    xyzview.setStyle({'sphere':{'radius':0.3},'stick':{'radius':0.1}})
    xyzview.setBackgroundColor('white')
    xyzview.addPropertyLabels("atom")
    xyzview.zoomTo()
    showmol(xyzview, height=500, width=800)

def optimize_geometry(atomic_numbers, positions):
    atoms = Atoms(atomic_numbers, 
                  positions)
    
    calc = LennardJones()
    atoms.calc = calc

    dyn = BFGS(atoms)
    dyn.run(fmax=0.05)

    updated_positions = atoms.get_positions()
    xyz_block = get_xyz_block(atomic_symbols = atoms.get_chemical_symbols(), 
                              positions = updated_positions)
    
    return xyz_block


def get_xyz_block(atomic_symbols: list[str], positions: np.ndarray) -> str:
    xyz_block = f"{len(atomic_symbols)}\n \n"
    for i, (symbol, position) in enumerate(zip(atomic_symbols , positions)):
        if i == len(atomic_symbols)-1:
            xyz_block += f"{symbol}    {position[0]:.6f}     {position[1]:.6f}     {position[2]:.6f}"
        else:
            xyz_block += f"{symbol}    {position[0]:.6f}     {position[1]:.6f}     {position[2]:.6f}\n"
    return xyz_block



atomic_numbers_to_symbols = {1: "H", 
                             2: "He", 
                             3: "Li", 
                             4: "Be", 
                             5: "B", 
                             6: "C", 
                             7: "N", 
                             8: "O", 
                             9: "F",
                             10: "Ne",
                             11: "Na", 
                             12: "Mg", 
                             13: "Al", 
                             14: "Si", 
                             15: "P", 
                             16: "S", 
                             17: "Cl", 
                             18: "Ar", 
                             19: "K",
                             20: "Ca",
                             21: "Sc", 
                             22: "Ti", 
                             23: "V", 
                             24: "Cr", 
                             25: "Mn", 
                             26: "Fe", 
                             27: "Co", 
                             28: "Ni", 
                             29: "Cu",
                             30: "Zn",
                             31: "Ga", 
                             32: "Ge", 
                             33: "As", 
                             34: "Se", 
                             35: "Br", 
                             36: "Kr", 
                             37: "Rb", 
                             38: "Sr", 
                             39: "Y",
                             40: "Zr",
                             41: "Nb", 
                             42: "Mo", 
                             43: "Tc", 
                             44: "Ru", 
                             45: "Rh", 
                             46: "Pd", 
                             47: "Ag", 
                             48: "Cd", 
                             49: "In",
                             50: "Sn",
                             51: "Sb", 
                             52: "Te", 
                             53: "I", 
                             54: "Xe", 
                             55: "Cs", 
                             56: "Ba", 
                             57: "La", 
                             58: "Ce", 
                             59: "Pr",
                             60: "Nd",
                             61: "Pm", 
                             62: "Sm", 
                             63: "Eu", 
                             64: "Gd", 
                             65: "Tb", 
                             66: "Dy", 
                             67: "Ho", 
                             68: "Er", 
                             69: "Tm",
                             70: "Yb",
                             71: "Lu", 
                             72: "Hf", 
                             73: "Ta", 
                             74: "W", 
                             75: "Re", 
                             76: "Os", 
                             77: "Ir", 
                             78: "Pt", 
                             79: "Au",
                             80: "Hg",
                             81: "Tl",
                             82: "Pb",
                             83: "Bi",
                             84: "Po",
                             85: "At"
                             }

