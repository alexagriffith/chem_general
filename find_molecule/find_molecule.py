import openbabel
import numpy as np
import re


def testopenbabel(reduced_xyz_filename,molecule):
    """
    Use open babel to convert the reduced xyz file ( file containing only atoms that are present in the input molecule )
    to a SMILES string, and outputs an smi string that matches the length of the input molecule and the xy coords that correspond to
    the smi string match.

    :param reduced_xyz_filename: xyz file containing symbol and x, y, z coords in a row --> H 0.0 0.0 1.0
    :param molecule: list of tuples representing the input molecule.. example: HSO3 --> [(H,1),(S,1),(O,3)]
    :return: molecule_str_smi: string matching the length of the input molecule: example HSOOO matches the tuple above
    :return: smi_xy_coords: X,Y coordinates in a list of tuples, matching molecule_str_smi  [(x,y),(x2,y2)]

    """

    # convert xyz to smiles
    conv = openbabel.OBConversion()
    conv.SetInAndOutFormats("xyz", "smi")
    mol = openbabel.OBMol()
    conv.ReadFile(mol, reduced_xyz_filename)
    #explicit hydrogens
    conv.AddOption("h", conv.OUTOPTIONS)
    # print x,y coords ( does not print any hydrogen coords if in bracket with another atom
    # example : [OH] hydrogen coord will not print, [H] hydrogen coord will print
    conv.AddOption("x", conv.OUTOPTIONS)
    conv.Convert()
    smi_data = (conv.WriteString(mol)).split()
    # split connectivity data by period ".", where each period separates connected atoms
    smi_string = (smi_data[0]).split(".")
    #print(smi_string)


    # get the number of atoms in the input molecule
    molecule_len = 0
    for atom in molecule:
        molecule_len += int(atom[1])

    # reduce the list of possible atoms by using smi string
    # find the smi string representations out of the smi output
    for index, set in enumerate(smi_string):
        if len(set) > molecule_len:
            result = "".join(i for i in re.findall(r"(?i)\b[a-z]+\b", set) if i.isalpha())
            if len(result) == molecule_len:
                # string that matches the input
               molecule_str_smi = result
                # index assuming everything in front of the molecule has one element example: [O].[S].[inputmolecule] index = 2
               mol_idx = index

    # if H in molecule_str_smi, update length for the indexing of the smi
    if 'H' in molecule_str_smi:
        # split on the H
        noH_molecule_stri_smi = molecule_str_smi.replace("H","")
        # assuming the hydrogens are always the last characters in the brackets
        noH_molstrsmi_len = len(noH_molecule_stri_smi)
        # will give correct range of molecule index, example: .[S][O][O][OH] has an XY coord in the range of 4, not 5 since H wont be in the list of coords



    count = 0
    # find the index in the smi string, starting from the beginning until the index that supposes everything before it is one element
    for i in smi_string[:mol_idx]:
        # if length of i in smi_string is greater than 3 ( two brackets + one characer, [O]), need to update to longer index
        if len(i) > 3:
            # if longer than 3 AND an H is present, the H will not count towards the index (no XY coord)
            if 'H' in i:
                #find_len = (re.findall("\[([^[\]]*)\]", i))
                # assumes H is always at the end
                no_H = i.split("H")[0]
                # len - 1 is accounting for the left bracket.. example: [OH] --> [O
                count += (len(no_H) - 1)
            else:
                # if no H, need to update index for every other element there is.
                # len - 2 where 2 is for the two brackets [SO] .. number of elements  = 4 - 2 = 2
                count += (len(i) - 2)
        else:
            count += 1
        # index range of molecule in smi string to get XY coordinats from indx1 to indx2
        smi_molecule_range = [count, (count + noH_molstrsmi_len)]



    # get xy coordinates in the smi_molecule_range
    xy_data = smi_data[2].split(",")
    xy_coords = []
    # from list [x1,y1,x2,y2] --> [(x1,y1),(x2,y2)]
    for idx, x in enumerate(xy_data):
        if idx % 2 == 0:
            # make x,y list
            xy_coord = (xy_data[idx],xy_data[idx+1])
            xy_coords.append(xy_coord)
        else:
            continue

    # from x,y list, get the coords within the index range of the molecule
    smi_xy_coords = []
    for xy_coord in xy_coords[smi_molecule_range[0]:smi_molecule_range[1]]:
       smi_xy_coords.append(xy_coord)


    return(molecule_str_smi,smi_xy_coords)


def Molecule(formula=None):
    """
    A simple molecule representation.

    Parameters
    __________

    formula: string input of the molecule wanted. examples: HSO3, H20, NO2.

    Returns
    _______

    mlelements: "molecule elements" a list of tuples giving information about the atom and number of times
    it's in the molecule  --> H2O = [(H,2),(O,1)]

    """
    mlelements = []

    if formula is None:
        formula = ''

    for index in range(len(formula)):
        char = formula[index]
        # if char in list is a character
        if char.isalpha():
            # if not at the end of the string,
            if (index + 1) != len(formula) and formula[index+1].isalpha() == False:
                mlelements.append((char, str(formula[index+1])))

            else:
                # if element X shows up only once  (X, 1)
                mlelements.append((char, str(1)))

        else:
            continue

    return (mlelements)


def check_unique_rep(mlelements=None):
    """
    check if the list of tuples is unique. If not, create a new list of tuples that is.
    Example [(H,1),(H,1),(O,2)] --> [(H,2),(O,1)]
    :param mlelements: string, list of tuples
    :return: unique_mlelements: string, list of tuples
    """

    if mlelements is None:
        pass

    unique_mlelements = []
    for x in set(mlelements):
        # count how many instances of each element
        counter = x,mlelements.count(x)
        if counter[1] != 1:
            # if more than one instance, create a new tuple with (element, number of instances)
            new_rep = (counter[0][0], str(counter[1]))
            unique_mlelements.append(new_rep)
        else:
            unique_mlelements.append(x)

    return(unqiue_mlelements)

def open_xyzfile(filename):
    """
    Open an xyz file, read the lines, and create an array

    :param filename: string, name of xyz file
    :return: string list, (id, symbol, (x,y,z)) example: (0, C, (1.0,1.0,0.0))
    """
    atoms_info = []
    f = open(filename, "r")
    total = f.readline()
    comment = f.readline()
    id = 0
    for line in f:
        line = line.split()
        atom_info = (id, line[0], (line[1:4]))
        atoms_info.append(atom_info)
        id += 1
    f.close()
    return(atoms_info)


def get_molecule(atoms_info=None, molecule=None):
    """
    :param atoms_info: list
    :param molecule: list
    :return: symbol_match, list of all elements in xyz file that match with the elements in the molecular formula
    """
    if atoms_info is None:
        atoms_info = []
    if molecule is None:
        molecule = Molecule()

    # reduce the xyz file (possible molecule coordinates) by searching the set of elements found in input molecule
    symbol_match = []
    for atom_symbol in atoms_info:
        for mol_symbol in molecule:
            if atom_symbol[1] == mol_symbol[0]:
                symbol_match.append(atom_symbol)
    return(symbol_match)




def calculate_distance(atom1,atom2):
    """
    calculate the distance btwn two atoms. This is a doc string
    explain what the function does but not how it does it. Dashes have to be exactly as long as the word.
    Parameters
    ----------
    atom1: list
        longer explanation of atom1. A list of coordinates [x, y, z].
    atom2: list
        list of coordinates [x, y, z].
    Returns
    -------
    bond_lengths: float
        The distance between atoms.
    Examples
    --------
    >> calculate_distance([0,0,0],[0,0,1])
    1.0
    """

    x_dist = float(atom1[0]) - float(atom2[0])
    y_dist = float(atom1[1]) - float(atom2[1])
    z_dist = float(atom1[2]) - float(atom2[2])
    distance = np.sqrt(x_dist**2 + y_dist**2 + z_dist**2)
    return (distance)


def bond_check(bond_distance, minimum_value=0.0, maximum_value=1.8): # min and max values = x, where x is the default value if not specified in the cell below.
    """
    This function checks if the distance between two atoms is within a min and max distance.
    Also we want to check if there is a real bond or not because the distance might not work sometimes.
    Parameters
    ----------
    bond_distance:  float
        The distance between two atoms
    minimum_value: float
        The minimum distance between two atoms
    maximum_value:
        The max distance between two atoms
    Returns
    -------
    True if bond
    False if not a bond
    """

    if bond_distance > minimum_value and bond_distance < maximum_value:
        return True
    else:
        return False

def match_xyz(symbol_match,smi_xy_coords):
    match_coords = []
    for coordsym in symbol_match:
        for coordsmi in smi_xy_coords:
            if str(round(float(coordsym[2][0]),4)) in coordsmi[0]:
                #print(coordsym,coordsmi)
                match_coords.append(coordsym)
    return(match_coords)



def match_H_xyz(symbol_match=None):
    """
    Isolate the hydrogen coordinates if the molecule has a hydrogen in it.
    This is necessary since x,y coordinates of hydrogens that are part of a molecule do not show up in a SMILES
    :param symbol_match: string, list of tuples containing elements that match with input molecule [id, symbol, (x,y,z)]
    :return: hydrogen_coords: string, list of tuples containing only hydrogens [id, symbol, (x,y,z)]
    """
    if symbol_match is None:
        pass

    # make an array with only hydrogens
    hydrogen_coords = []
    for x in symbol_match:
        if x[1] == 'H':
            hydrogen_coords.append(x)
    return(hydrogen_coords)


def make_xyz(filename=None,coords):
    """
    Makes an xyz file from a list of tuples
    :param coords: string, list of tuples
    :return: makes file
    """

    if filename is None:
        pass

    with open(filename, "w") as f:
        f.write(str(len(coords)) + '\n')
        f.write('\n')
        for id, symbol, [x, y, z] in coords:
            f.write(f"{symbol} {x} {y} {z}\n")


def make_exclusive_array(atoms_info,match_coords):
    """
    Takes atoms_info, which is a list [id, symbol, (x,y,z)] of all the atoms from the original xyz file
    and makes a list of everything NOT belonging to the molecule.
    For interaction energies.

    :param atoms_info: string, list of tuples
    :param match_coords: string, list of tuples
    :return:  exclusive_atoms: string, list of tuples

    """
    # make an array with everything BUT the coords in match_coords
    exclusive_atoms = []
    for x in atoms_info:
        if x not in match_coords:
            exclusive_atoms.append(x)
    return (exclusive_atoms)


if __name__ == "__main__":
    """
    The main goal of this script is to automatically find a molecule within an xyz file 
    and make two xyz coords: one with the molecule and the other containing everything but the molecule. 
    
    Can be used for the geometries needed to calculate interaction energies, etc. 
    Example: get xyz of a substrate in an active site 
    
    REQUIREMENTS: (see requirements.txt) mainly need open babel installed and python 3 
    
    """
    # atoms_info is a list of (index, symbol, (x,y,z))
    input_molecule = 'HSO3'
    xyz_filename = "testlarge.xyz"
    atoms_info = open_xyzfile("testlarge.xyz")
    # molecule is a list of tuples
    molecule = Molecule(input_molecule)
    symbol_match = get_molecule(atoms_info=atoms_info, molecule=molecule)

    # write an xyz file with reduced list of atoms that match with the input molecule
    make_xyz("symbol_match.xyz",symbol_match)

    # find a match from smi string representation of the bound atoms in the given xyz file, example: H2O --> HHO
    # smi_xy_coords are sets of X,Y coordinates corresponding to the match_from_smi string
    match_from_smi, smi_xy_coords = testopenbabel("symbol_match.xyz",molecule=molecule)

    # get a tuple list of molecule example: H2O  [(H,1),(H,1),(O,1)]
    extracted_smi_molecule = Molecule(str(match_from_smi))
    # get a unique representation of the smi string example: H2O [(H,2),(O,1)]
    unique_smi_rep = check_unique_rep(extracted_smi_molecule)

    # test if the input and the smi unique tuple list match
    if set(molecule) == set(unique_smi_rep):
        # find the matching xyz coords
        match_coords = match_xyz(symbol_match,smi_xy_coords)
        count = 0
        for coordA in match_coords:
            for coordB in match_coords:
                if coordA[0] != coordB[0]:
                    distance = calculate_distance(coordA[2],coordB[2])
                    if bond_check(distance):
                        count += 1
                        if count == len(match_coords):
                            pass
        print("Y")
    else:
        print("Your molecule is lonely! :( No match was found.")

    # find hydrogen coords
    hydrogen_coords = match_H_xyz(symbol_match)

    # find hydrogen(s) in the molecule
    for coordC in match_coords:
        for coordD in hydrogen_coords:
            distance = calculate_distance(coordC[2],coordD[2])
            if bond_check(distance,maximum_value=1.1):
                match_coords.append(coordD)

    # make two xyz files, one containing input molecule and the other containing everything else
    make_xyz(f"{input_molecule}.xyz",match_coords)
    exclusive_atoms = make_exclusive_coords(atoms_info,match_coords)
    make_xyz("exclusive_atoms.xyz",exclusive_atoms)



