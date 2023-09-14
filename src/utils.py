import warnings
from typing import Union

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
    return {"n_mols": int, "n_atoms_min": int, "n_atoms_max": int, "atomic_numbers": list, "unit_cell": list, "min_radius": float, "max_radius": float, "random_seed": int, "name_id": str}