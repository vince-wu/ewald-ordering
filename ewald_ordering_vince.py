from pymatgen.alchemy.transmuters import StandardTransmuter
from pymatgen.transformations.standard_transformations import OrderDisorderedStructureTransformation, \
    OxidationStateDecorationTransformation
from pymatgen.alchemy.materials import TransformedStructure
from pymatgen.io.cif import CifWriter
from pymatgen.io.vasp import Poscar
from pymatgen import Structure
import tkinter as tk
from tkinter import filedialog
import os

#Prompt File Selection
root = tk.Tk()
root.withdraw()
path = filedialog.askopenfilename()

#Read cif file
cryst = TransformedStructure(Structure.from_file(path))

#Create Oxidation State Transform - change this line to assign oxidation states!
oxidation_transform = OxidationStateDecorationTransformation({"Fe": 3, "Mn": 2, "O": -2, "V": 5, "Na": 1})

#Create Ewald Ordering Transform
ordering_transform = OrderDisorderedStructureTransformation()

#Run transformations on disordered crystal structure
transmuter = StandardTransmuter([cryst],[oxidation_transform, ordering_transform])
newCryst = transmuter.transformed_structures[0].final_structure

print("Ewald optimization successful!")

#Save to CIF
filename = os.path.splitext(os.path.basename(path))[0] + '_ewald'
w = CifWriter(newCryst)
w.write_file(filename + '.cif')
print("Cif file saved to %s.cif"%(filename))

#Save to POSCAR
poscar = Poscar(newCryst)
poscar.write_file(filename)
print("POSCAR file saved to %s"%(filename))

