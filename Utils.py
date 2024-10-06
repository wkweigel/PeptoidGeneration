import random
from Reactions import Reaction

def generate_peptoid_strings(n:int, length:int, monomers:list):
    '''Generate a set of n random strings using a list of monomers
    Parameters:
        n: An integer representing the number of strings to generate .
        length: The length of the desired peptoids.
        monomers: A list of single-character monomers. 

    Returns:
        unique_strings: A list of peptoid strings of length n.
    '''
    
    unique_strings = set() #Using a set so there are no repeat strings
    
    while len(unique_strings) < n:
        random_string = ''.join(random.choices(monomers, k=length))
        unique_strings.add(random_string)
    
    return list(unique_strings)

def peptoid_smiles_from_string(peptoid_string:str, peptoid_dict:dict):
    '''Convert peptoid string into its corresponding smiles.
    Parameters:
        peptoid_string: An string representing the monomers in the desired peptoid. 
        peptoid_dict: A dictionary of corresponding monomer codes and smiles.
        
    Returns:
        final_peptoid: The smiles string of the fully assembled peptoid.
    '''
    #Start with the specified brominated monomer 'Z'
    #peptoid=Synthesis(peptoid_dict['Z'])

    #Start with protected bromoacetic acid
    peptoid=Reaction(int_halide=peptoid_dict['B'])

    for idx, monomer in enumerate(peptoid_string): #For cases involving
        if monomer == 'X':
            if idx==0:
                peptoid.int_halide_smiles=peptoid.amine_displacement(peptoid_dict[monomer],peptoid.int_halide_smiles)
            else:
                peptoid.int_halide_smiles=peptoid.amine_displacement(peptoid_dict[monomer],peptoid.int_halide_smiles)   
        else:
            if idx==0:
                peptoid.int_amine_smiles=peptoid.amine_displacement(peptoid_dict[monomer],peptoid.int_halide_smiles)
                peptoid.int_halide_smiles=peptoid.amide_bond_formation2(peptoid.int_amine_smiles)
            else:
                peptoid.int_amine_smiles=peptoid.amine_displacement(peptoid_dict[monomer],peptoid.int_halide_smiles)
                peptoid.int_halide_smiles=peptoid.amide_bond_formation2(peptoid.int_amine_smiles)


    if peptoid_string[-1] == 'X':
        final_peptoid=peptoid.deprotect_peptoid(peptoid_string, peptoid.int_halide_smiles)
    else:
        final_peptoid=peptoid.deprotect_peptoid(peptoid_string, peptoid.int_amine_smiles)

    return(final_peptoid)