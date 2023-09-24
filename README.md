# Generate XYZ files for molecules with random atomic coordinates
This app allows you to generate random structures with specified number of atoms, elements, unit cell size, and interatomic distances.
The created molecules can be optimized using GFN1-xTB.  
## Demo
<p align="center">
  <img src="https://github.com/ineporozhnii/mol_gen_xyz/blob/main/assets/mol_gen_demo.gif" alt="animated" />
</p>

## Run locally
### Create environment
 - `conda mol_gen_env create -f environment.yml`
 - `conda activate mol_gen_env`
### Run app
- `streamlit run app.py`
