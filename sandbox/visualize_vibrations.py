import streamlit as st
import py3Dmol
import numpy as np

# Function to read and parse the g98.out file
def read_g98_file(uploaded_file):
    coordinates = []
    frequencies = []
    displacements = []
    reading_coordinates = False
    reading_frequencies = False

    content = uploaded_file.read().decode('utf-8')  # Read and decode the file content
    lines = content.splitlines()
    lines_iter = iter(lines)  # Create an iterator from the list

    for line in lines_iter:
        if 'Standard orientation:' in line:
            reading_coordinates = True
            coordinates = []  # Reset coordinates to capture the latest block
            # Skip the next four lines (headers)
            for _ in range(4):
                next(lines_iter, None)
            continue

        if reading_coordinates:
            if '----' in line:
                reading_coordinates = False
                continue
            parts = line.split()
            if len(parts) == 6 and parts[0].isdigit():  # Ensure line has 6 parts and starts with a number
                coordinates.append({
                    'element': int(parts[1]),
                    'x': float(parts[3]),
                    'y': float(parts[4]),
                    'z': float(parts[5])
                })

        if 'Frequencies --' in line:
            reading_frequencies = True
            freqs = line.split()[2:]
            frequencies.extend([float(freq) for freq in freqs])

        if reading_frequencies and 'Atom AN' in line:
            reading_frequencies = False
            # Read the displacement vectors
            mode_disp = []
            for i in range(len(coordinates)):
                disp_line = next(lines_iter).split()[2:]
                mode_disp.append([float(d) for d in disp_line])
            displacements.append(mode_disp)

    return coordinates, frequencies, displacements

# Function to visualize the molecular structure using py3Dmol
def visualize_molecule(coordinates):
    mol_view = py3Dmol.view(width=800, height=400)
    xyz_str = f"{len(coordinates)}\n\n"
    for atom in coordinates:
        element = atom['element']
        x, y, z = atom['x'], atom['y'], atom['z']
        xyz_str += f"{element} {x:.4f} {y:.4f} {z:.4f}\n"
    
    mol_view.addModel(xyz_str, "xyz")
    mol_view.setStyle({"stick": {}})
    mol_view.zoomTo()
    return mol_view

# Function to visualize vibrational modes
def visualize_vibrations(coordinates, displacements, frequency_index):
    mol_view = py3Dmol.view(width=800, height=400)
    xyz_str = f"{len(coordinates)}\n\n"
    for atom in coordinates:
        element = atom['element']
        x, y, z = atom['x'], atom['y'], atom['z']
        xyz_str += f"{element} {x:.4f} {y:.4f} {z:.4f}\n"
    
    mol_view.addModel(xyz_str, "xyz")
    mol_view.setStyle({"stick": {}})
    mol_view.zoomTo()
    
    # Visualize vibration
    displacement = np.array(displacements[frequency_index]).reshape(-1, 3)
    for step in range(-10, 11, 2):
        scale = step * 0.1  # Scale for displacement
        vib_str = f"{len(coordinates)}\n\n"
        for i, atom in enumerate(coordinates):
            x = atom['x'] + scale * displacement[i][0]
            y = atom['y'] + scale * displacement[i][1]
            z = atom['z'] + scale * displacement[i][2]
            element = atom['element']
            vib_str += f"{element} {x:.4f} {y:.4f} {z:.4f}\n"
        mol_view.addModel(vib_str, "xyz")
#        mol_view.setStyle({"stick": {}})
        mol_view.setStyle({'stick': {}, 'sphere': {'radius': 0.5}})
    
    return mol_view

# Streamlit interface
st.title('Molecular Vibrations and Structure Visualization')

uploaded_file = st.file_uploader("Upload g98.out file", type=['out', 'txt'])

if uploaded_file is not None:
    # Read the file and parse data
    coordinates, frequencies, displacements = read_g98_file(uploaded_file)

    # Visualize the molecule
    st.subheader("Molecular Structure")
    mol_view = visualize_molecule(coordinates)
    st.components.v1.html(mol_view._make_html(), width=800, height=400)

    # Show frequencies and allow user to select one for visualization
    st.subheader("Vibrational Frequencies")
    if frequencies:
        st.write(f"Available Frequencies (cm^-1): {frequencies}")
        frequency_index = st.selectbox("Select a frequency to visualize the vibration:", range(len(frequencies)), format_func=lambda x: f"{frequencies[x]:.2f} cm^-1")
        vib_view = visualize_vibrations(coordinates, displacements, frequency_index)
        st.components.v1.html(vib_view._make_html(), width=800, height=400)
    else:
        st.write("No frequencies found in the file.")
else:
    st.write("Please upload a g98.out file to visualize.")

