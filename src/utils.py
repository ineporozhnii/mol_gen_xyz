import warnings
import py3Dmol
from rdkit import Chem
from typing import Union
from stmol import showmol
from ase import Atoms
from ase.optimize import BFGS
from tblite.ase import TBLite
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
    xyzview.setBackgroundColor('white')#('0xeeeeee')
   # xyzview.add_label(label_type="atomname")
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

    xyz_block = f"{len(atoms.get_atomic_numbers())}\n \n"
    for i, (symbol, position) in enumerate(zip(atoms.get_chemical_symbols() , updated_positions)):
        if i == len(atoms.get_atomic_numbers())-1:
            xyz_block += f"{symbol}    {position[0]:.5f}     {position[1]:.5f}     {position[2]:.5f}"
        else:
            xyz_block += f"{symbol}    {position[0]:.5f}     {position[1]:.5f}     {position[2]:.5f}\n"
    return xyz_block


atomic_numbers_to_symbols = {1: "H", 
                             2: "He", 
                             3: "Li", 
                             4: "Be", 
                             5: "B", 
                             6: "C", 
                             7: "N", 
                             8: "O", 
                             9: "F"}

