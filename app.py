import streamlit as st
from src.main import main

st.set_page_config(layout="wide", page_title="Generate random molecules")


def get_input():
    from src.utils import atomic_numbers_to_symbols
    atomic_symbols_to_numbers = dict((atomic_numbers_to_symbols[key], key) for (key, value) in atomic_numbers_to_symbols.items())

    n_mols = st.sidebar.slider("Select number of molecules to generate", 
                               min_value=1, 
                               max_value=100)
    
    
    n_atoms_min, n_atoms_max = st.sidebar.slider("Select minimum and maximum number of atoms in a molecule", 
                                                 min_value=2, 
                                                 max_value=20, 
                                                 value=(5, 10))
    
    
    atomic_symbols = st.sidebar.multiselect("Select allowed elements", list(atomic_numbers_to_symbols.values()), ["H", "C", "O", "F"])
    atomic_numbers = [atomic_symbols_to_numbers[symbol] for symbol in atomic_symbols]

    # Unit cell 
    unit_cell_x = st.sidebar.number_input('Unit cell size x axis', value=5.0)
    unit_cell_y = st.sidebar.number_input('Unit cell size y axis', value=5.0)
    unit_cell_z = st.sidebar.number_input('Unit cell size z axis', value=5.0)
    unit_cell = [unit_cell_x, unit_cell_y, unit_cell_z]

    min_distance, max_distance = st.sidebar.slider("Select minimum and maximum distance between atoms", 
                                                 min_value=0.1, 
                                                 max_value=30.0, 
                                                 value=(2.0, 10.0))

    random_seed = st.sidebar.number_input("Select random seed", 42)

    no_disconnected_mols = st.sidebar.checkbox('No disconnected molecules', value=True)

    return n_mols, n_atoms_min, n_atoms_max, atomic_symbols, atomic_numbers, unit_cell, min_distance, max_distance, random_seed, no_disconnected_mols


def app():

    n_mols, n_atoms_min, n_atoms_max, atomic_symbols, atomic_numbers, unit_cell, min_distance, max_distance, random_seed, no_disconnected_mols = get_input()

    generate_button = st.sidebar.button("Generate", type='primary', use_container_width=True)

    mol_list = []
    if generate_button:
        with st.spinner("Generating molecules..."):
            mol_list = main(n_mols=n_mols,
                n_atoms_min = n_atoms_min, 
                n_atoms_max = n_atoms_max,
                atomic_numbers = atomic_numbers, 
                unit_cell = unit_cell,
                min_distance = min_distance,
                max_distance = max_distance,
                random_seed = random_seed,
                no_disconnected_mols = no_disconnected_mols,
                return_mols = True)

            mol_dict = dict((f'generated_mol_{i}', mol) for (i, mol) in enumerate(mol_list))
            st.session_state["mol_dict"] = mol_dict

    if "mol_dict" in st.session_state:
        n_generated_molecules = len(st.session_state["mol_dict"])
        st.title(f"Generated {n_generated_molecules} molecules")
        selected_mol_id = st.selectbox("Selecte generated molecule", list(st.session_state["mol_dict"].keys()))

        col1, col2 = st.columns(2)
        with col1:
            st.header("Generated molecule:")
            st.write("Plot geometry")

        with col2:
            st.header("Geometry xyz:")
            st.text(st.session_state["mol_dict"][selected_mol_id].xyz_block)





app()
