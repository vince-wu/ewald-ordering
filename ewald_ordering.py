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


BASE_DIR = 'c:/Users/Vincent/Box/Clement Research/'

def generate_ewald_orderings(path, choose_file, oxidation_states, num_structures):
    """
    DESCRIPTION: Given a disordered CIF structure with at least one crystallographic site that is shared by more than one element, all permutations
                 will have their electrostatic energy calculated via an Ewald summation, given that all ion charges are specified. Ordered CIF structures will 
                 be generated, postpended with a number that indicates the stability ranking of the structure. For example, if a CIF file called 
                 "Na2Mn2Fe(VO4)3.cif" is inputted with num_structures=3, then the function will generate 3 ordered output files, 
                 "Na2Mn2Fe(VO4)3-ewald-1", "Na2Mn2Fe(VO4)3-ewald-2", and "Na2Mn2Fe(VO4)3-ewald-3", with "Na2Mn2Fe(VO4)3-ewald-1" being the most stable and
                 "Na2Mn2Fe(VO4)3-ewald-3" being the least stable. 
    PARAMETERS:
        path: string
            The file path to the CIF file. The CIF file must be a disordered structure, or else an error will occur. A disordered structure will have one
            site that is occupied by more than one element, with occupancies less than 1: i.e. Fe at 0,0,0.5 with an occupancy of 0.75, and Mn at the same 
            site 0,0,0.5 with an occupancy of 0.25. This function cannot handle multivalent elements, for example it cannot handle a structure that has Mn
            in both the 2+ and 3+ redox state.
        choose_file: boolean
            Setting this parameter to True brings up the file explorer dialog for the user to manually select the CIF file
        oxidation_states: dictionary
            A dictionary that maps each element in the structure to a particular oxidation state. E.g. {"Fe": 3, "Mn": 2, "O": -2, "V": 5, "Na": 1, "Al":3}
            Make sure that all elements in the structure is assigned an oxidation state. It is ok to add more elements than needed.
        num_structures: int
            The number of strcutures to be outputted by the function. There can be hundreds or thousands of candidate structures, however in practice only
            the first few (or even only the first, most stable) structures are needed.
    RETURNS: None
    """
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

    #Create Oxidation State Transform 
    oxidation_transform = OxidationStateDecorationTransformation(oxidation_states)

    #Create Ewald Ordering Transform object which will be passed into the transmuter
    ordering_transform = OrderDisorderedStructureTransformation(symmetrized_structures=True)

    # apply the order-disorder transform on the structure for any site that has fractional occupancies
    transmuter = StandardTransmuter([cryst],[oxidation_transform, ordering_transform], extend_collection=num_structures)
    print("Ewald optimization successful!")
    for i in range(num_structures):
        newCryst = transmuter.transformed_structures[i].final_structure
        #Save to CIF
        structure_name = os.path.splitext(os.path.basename(path))[0]
        save_directory = structure_name
        if not os.path.isdir(save_directory):
            os.mkdir(save_directory)
        filename = structure_name + '/' + structure_name + '-ewald' + '-%i'%(i+1)
        w = CifWriter(newCryst)
        w.write_file(filename + '.cif')
        print("Cif file saved to {}.cif".format(filename))
        #Save to POSCAR
        poscar = Poscar(newCryst)
        poscar.write_file(filename)
        print("POSCAR file saved to {}".format(filename))

generate_ewald_orderings(
    path = BASE_DIR + 'Python Scripts/ewald-ordering/Na2Mn2Fe(VO4)3-disordered.cif',
    choose_file=False,
    oxidation_states={"Fe": 3, "Mn": 2, "O": -2, "V": 5, "Na": 1, "Al":3, "P":5},
    num_structures=16
)



