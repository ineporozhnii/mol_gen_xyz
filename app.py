import streamlit as st
from src.main import main
from src.utils import get_connected_mol, display_molecule, optimize_geometry


st.set_page_config(layout="wide", page_title="Generate random molecules")


def get_input():
    st.sidebar.title("Molecule Generator")
    from src.utils import atomic_numbers_to_symbols
    atomic_symbols_to_numbers = dict((atomic_numbers_to_symbols[key], key) for (key, value) in atomic_numbers_to_symbols.items())

    n_mols = st.sidebar.slider("Select number of molecules to generate", 
                               min_value=1, 
                               max_value=100,
                               value=5)
    
    
    n_atoms_min, n_atoms_max = st.sidebar.slider("Select minimum and maximum number of atoms in a molecule", 
                                                 min_value=2, 
                                                 max_value=20, 
                                                 value=(2, 5))
    
    
    atomic_symbols = st.sidebar.multiselect("Select allowed elements", 
                                            [f"{str(number)}: {symbol}" for (number, symbol) in atomic_numbers_to_symbols.items()], 
                                            ["1: H", "6: C", "8: O", "9: F"])
    atomic_numbers = [int(symbol.split(":")[0]) for symbol in atomic_symbols]

    # Unit cell 
    unit_cell_x = st.sidebar.number_input('Unit cell size x axis', value=3.0)
    unit_cell_y = st.sidebar.number_input('Unit cell size y axis', value=3.0)
    unit_cell_z = st.sidebar.number_input('Unit cell size z axis', value=3.0)
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
    #st.sidebar.write("**Created by Ihor Neporozhnii**, [GitHub](https://github.com/ineporozhnii/mol_gen_xyz)")
    st.sidebar.markdown("<h3 style='text-align: center; color: #31333F;'>Created by Ihor Neporozhnii</h3>", unsafe_allow_html=True)
    st.sidebar.markdown("""<a style='display: block; text-align: center;' href="https://github.com/ineporozhnii/mol_gen_xyz">GitHub</a>""", unsafe_allow_html=True)

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
            
            if mol_list is None:
                st.error("Generation failed, please try with different parameters")
            else:
                mol_dict = dict((f'generated_mol_{i}', mol) for (i, mol) in enumerate(mol_list))
                st.session_state["mol_dict"] = mol_dict

    if "mol_dict" in st.session_state:
        n_generated_molecules = len(st.session_state["mol_dict"])
        if n_generated_molecules == 1:
            st.title(f"Generated {n_generated_molecules} molecule")
        else:
            st.title(f"Generated {n_generated_molecules} molecules")
        selected_mol_id = st.selectbox("Selecte generated molecule", list(st.session_state["mol_dict"].keys()))

        with st.expander('Generated molecule', expanded=True):
            col1, col2 = st.columns(2)
            with col1:
                st.header("Generated molecule:")
                mol_block = get_connected_mol(st.session_state["mol_dict"][selected_mol_id].xyz_block)
                display_molecule(mol_block)
            with col2:
                st.header("Generated xyz:")
                st.text(st.session_state["mol_dict"][selected_mol_id].xyz_block)
                st.download_button('Download geometry XYZ file', 
                                    st.session_state["mol_dict"][selected_mol_id].xyz_block, 
                                    file_name=f'{selected_mol_id}_geometry.xyz')

        with st.expander("Optimized molecule", expanded=True):
            optimize_button = st.button("Pre-optimize geometry with BFGS + Lennard-Jones", type='primary', use_container_width=True)
            col1, col2 = st.columns(2)
            with col1:
                if optimize_button:
                    optimized_xyz_block = optimize_geometry(st.session_state["mol_dict"][selected_mol_id].atomic_numbers, 
                                                        st.session_state["mol_dict"][selected_mol_id].positions)
                
                    mol_block = get_connected_mol(optimized_xyz_block)

                    st.header("Optimized molecule:")
                    display_molecule(mol_block)
            with col2:
                if optimize_button:
                    st.header("Optimized xyz:")
                    st.text(optimized_xyz_block)
                    st.download_button('Download pre-optimized geometry XYZ file', 
                                       optimized_xyz_block, 
                                       file_name=f'{selected_mol_id}_geometry_pre_optimized.xyz')

        st.divider()
        save_all_button = st.button("Save all generated XYZs", type='primary', use_container_width=True)
        if save_all_button:
            st.write("saving")





app()
