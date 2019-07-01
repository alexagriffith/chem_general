# chem_general
general scripts useful for quantum chemistry programs


find_molecule/
 The main goal of this script is to automatically find a molecule within an xyz file 
    and make two xyz coords: one with the molecule and the other containing everything but the molecule. 
    
    Can be used for the geometries needed to calculate interaction energies, etc. 
    Example: get xyz of a substrate in an active site 
    
    REQUIREMENTS: (see requirements.txt) mainly need open babel installed and python 3 
    
