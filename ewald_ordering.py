from pymatgen.alchemy.transmuters import StandardTransmuter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
from pymatgen.transformations.standard_transformations import OrderDisorderedStructureTransformation, \
    OxidationStateDecorationTransformation
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp import Poscar
from pymatgen.core import Structure
import tkinter as tk
from tkinter import filedialog
import os


BASE_DIR = 'c:/Users/Vincent Wu/Box/Clement Research/' # the directory all .res files are in
if not os.path.exists(BASE_DIR):
    BASE_DIR = 'c:/Users/Vincent/Box/Clement Research/' # alternate directory for a different device
    
OUTPUT_DIR = BASE_DIR + 'Python Scripts/ewald-ordering/'


def generate_ewald_orderings(disordered_cif_path, choose_file, oxidation_states, num_structures):
    # open file dialog if file is to be manually chosen
    if choose_file:
        root = tk.Tk()
        root.withdraw()
        path = filedialog.askopenfilename()
    #Read cif file
    cryst = Structure.from_file(path)
    analyzer = SpacegroupAnalyzer(cryst)
    space_group = analyzer.get_space_group_symbol()
    print(space_group)
    symm_struct = analyzer.get_symmetrized_structure()
    cryst = TransformedStructure(symm_struct)

    #Create Oxidation State Transform - change this line to assign oxidation states!
    oxidation_transform = OxidationStateDecorationTransformation({"Fe": 3, "Mn": 2, "O": -2, "V": 5, "Na": 1
    , "Al":3, "P":5})


    #Create Ewald Ordering Transform
    ordering_transform = OrderDisorderedStructureTransformation(symmetrized_structures=True)

    #Run transformations on disordered crystal structure
    num_structures = 2 #Specify the number of structures to output
    transmuter = StandardTransmuter([cryst],[oxidation_transform, ordering_transform], extend_collection=num_structures)
    print("Ewald optimization successful!")
    for i in range(num_structures):
        newCryst = transmuter.transformed_structures[i].final_structure
        #Save to CIF
        filename = os.path.splitext(os.path.basename(path))[0] + '_ewald' + '_%i'%(i+1)
        w = CifWriter(newCryst)
        w.write_file(filename + '.cif')
        print("Cif file saved to %s.cif"%(filename))

        #Save to POSCAR
        poscar = Poscar(newCryst)
        poscar.write_file(filename)
        print("POSCAR file saved to %s"%(filename))




