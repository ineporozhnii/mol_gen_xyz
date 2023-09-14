import os
import json
import random
import warnings
import numpy as np
from utils import check_input
from generator import generate_mol

def main(n_mols: int = 10, 
         n_atoms_min: int = 5, 
         n_atoms_max: int = 10, 
         atomic_numbers: list[int] = [1, 6, 8], 
         unit_cell: list[float] = [5.0, 5.0, 5.0],
         min_radius: float = 2.0,
         max_radius: float = 10.0,
         random_seed: int = 42,
         name_id: str = "mol",
         generated_limit: int = 10
         ):
    
    random.seed(random_seed)

    mol_id_list = []
    all_atoms_list = []
    

    n_accepted_mols = 0
    n_generated_mols = 0

    while n_accepted_mols < n_mols:
    

        mol = generate_mol()

        n_generated_mols += 1
        if n_generated_mols == generated_limit and n_accepted_mols == 0:
            warnings.warn(f"Generated: {n_generated_mols}, accepted: 0, stopping.", stacklevel=3)
            break

    return

if __name__ == '__main__':
    if os.path.exists('config.json'):
        try: 
            with open('config.json') as config_file:
                config_dict = json.load(config_file)
        except Exception as e:
            print(e)
            print("! Error reading config.json. Please, make sure the config file is correct.")
            print("! Proceeding with default parameters.")
            config_dict = {}
    else:
        warnings.warn("No config.json file in src/ directory. \n! Proceeding with default parameters.", stacklevel=3)
        config_dict = {}
    
    clean_config_dict = check_input(config_dict=config_dict)
    main(**clean_config_dict)
