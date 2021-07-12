from pymatgen.alchemy.transmuters import StandardTransmuter
from pymatgen.core import structure
from pymatgen.core.structure import IStructure
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
                 "Na2Mn2Fe(VO4)3-ewald-3" being the least stable. Note that this function does not take into account the symmetry of the crystal, and thus
                 it may give several structures which are symmetrically identical under a space group. Use "find_unique_structures" to isolate unique orderings.
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
    num_structures = len(transmuter.transformed_structures)
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

def find_unique_structures(base_crystal_path, directory):
    """
    DESCRIPTION: Given a base crystal structure along with a directory containing CIF structures, determine how many unique configurations are in the directory,
                 and assign a configuration ID number to each CIF structure in the directory. Results can be seen in the outputted unique.csv file. The function
                 uses the base crystal structure provided to find the space group and associated symmetry operations, which are then used to see whether the CIF
                 structures in the directory are equivalent. 
                 
    PARAMETERS:
        base_crystal_path: string
            The path to the base crystal CIF structure. Note that this structure should contain the space group that you want to check all other structures to.
            If the symmetry of this structure is incorrect, you will get unexpected results. 
        directory: string
            The path to the directory containing the CIF structures to compare. Note that this function assumes that the CIF names are all prefixed by a dash
            and followed by an interger, ie "LiCoO2-1", "LiCoO2-2", "LiCoO2-3"..., ect.
    RETURNS: None
    """
    # convert all files (they should all be cif files) in the directory into Structures and place them into a list
    cif_file_names = [os.path.splitext(os.path.basename(path))[0] for path in os.listdir(directory) if path.endswith('.cif')]
    cif_files = [directory + path for path in os.listdir(directory) if path.endswith('.cif')]
    cif_files = sorted(cif_files, key = lambda x: int(x.split("-")[-1].replace(".cif", "")))
    cif_file_names = sorted(cif_file_names, key = lambda x: int(x.split("-")[-1]))
    structures = [IStructure.from_file(path) for path in cif_files]
    print(cif_file_names)
    # Get the base crystal structure
    base_crystal = IStructure.from_file(base_crystal_path)

    # Get the space group of the base crystal
    analyzer = SpacegroupAnalyzer(base_crystal)
    spacegroup = analyzer.get_space_group_operations()
    print(spacegroup)


    id = 1
    num_structures = len(structures)
    config_dict = {}
    unique_structs = []
    for i in range(num_structures):
        config1 = structures[i]
        unique_cif_name = cif_file_names[i]
        if not unique_cif_name in config_dict:
            config_dict[unique_cif_name] = id
            print("%s Unique Structures Found..."%(id))
            unique_structs.append(config1)
            id += 1
            for j in range(num_structures):
                cif_name = cif_file_names[j]
                config2 = structures[j]
                if not cif_name in config_dict:
                    isEquivalent = spacegroup.are_symmetrically_equivalent(config1, config2)
                    if isEquivalent:
                        config_dict[cif_name] = config_dict[unique_cif_name]
    print("----------------------------------------------")
    print("%s Total Unique Structures Found"%(id-1))
    print("----------------------------------------------")
    with open('unique.csv', 'w') as f:
        for key in config_dict.keys():
            f.write("%s,%s\n"%(key,config_dict[key]))
    crystal_name = os.path.splitext(os.path.basename(base_crystal_path))[0]
    directory_path = "{}-unique".format(crystal_name)
    if not os.path.isdir(directory_path):
        os.mkdir(directory_path)

    num_unique_structures = len(unique_structs)
    for i in range(num_unique_structures):
        config = unique_structs[i]
        #Save to CIF
        filename = crystal_name + '_ewald' + '_{}'.format(i+1)
        w = CifWriter(config)
        w.write_file(directory_path + '/' + filename + '.cif')
        # print("Cif file saved to {}.cif".format(filename))

        #Save to POSCAR
        poscar = Poscar(config)
        poscar.write_file(directory_path + '/' + filename)
        # print("POSCAR file saved to {}".format(filename))

base_disordered_crystal_path = BASE_DIR + 'Python Scripts/ewald-ordering/Na2Mn2Fe(VO4)3-disordered.cif'
# generate_ewald_orderings(
#     path = base_disordered_crystal_path,
#     choose_file=False,
#     oxidation_states={"Fe": 3, "Mn": 2, "O": -2, "V": 5, "Na": 1, "Al":3, "P":5},
#     num_structures=1000
# )

find_unique_structures(
    base_crystal_path=base_disordered_crystal_path,
    directory=BASE_DIR + "Python Scripts/ewald-ordering/Na2Mn2Fe(VO4)3-disordered/"
)

