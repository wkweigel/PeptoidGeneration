
from rdkit import Chem
from rdkit.Chem import AllChem

class Reaction:
    '''A Class for iteratively assembling peptoids from amine monomers'''
    def __init__(self, inital_dipeptide=None, int_halide=None ):
        #Set initial synthesis variables
        if int_halide==None:
            self.inital_monomer=inital_dipeptide
            self.int_halide_smiles=''

        if inital_dipeptide==None:
            self.inital_monomer=''
            self.int_halide_smiles=int_halide

        self.int_amine_smiles=''
        self.int_deprotected_smiles=''
        self.int_deprotected_mol=None

        #Define the reaction SMARTS
        self.displacement_smarts = '[N;!H0!$(NC=[!#6;!P]):1].[#6:2][F,Cl,Br,I]>>[N:1][#6:2]'
        self.amidation_smarts = '[C:1](=[O:2])-[OD1].[N!H0:3]>>[C:1](=[O:2])[N:3]'
        self.boc_deprotection_smarts='[N:1]C(=O)OC(C)(C)(C)>>[N:1]'
        self.gly_deprotection_smarts='[N:1]C(=P)OC(C)(C)(C)>>[N:1]'
        self.Ar_boc_deprotection_smarts='[n:1]C(=O)OC(C)(C)(C)>>[nH:1]'
        self.ester_deprotection_smarts='[CX3:1](=[OX1:2])-[OX2:3]-[C](C)(C)(C)>>[C:1](=[O:2])-[O:3]'
        self.alcohol_deprotection_smarts='[*:1]-[OX2:3]-[C](C)(C)(C)>>[*:1]-[O:3]'

        #Create reactions from the SMARTS
        self.displacement_reaction=AllChem.ReactionFromSmarts(self.displacement_smarts)
        self.amidation_reaction=AllChem.ReactionFromSmarts(self.amidation_smarts)
        self.boc_deprotection_reaction=AllChem.ReactionFromSmarts(self.boc_deprotection_smarts)
        self.gly_deprotection_reaction=AllChem.ReactionFromSmarts(self.gly_deprotection_smarts)
        self.Ar_boc_deprotection_reaction=AllChem.ReactionFromSmarts(self.Ar_boc_deprotection_smarts)
        self.ester_deprotection_reaction=AllChem.ReactionFromSmarts(self.ester_deprotection_smarts)
        self.alcohol_deprotection_reaction=AllChem.ReactionFromSmarts(self.alcohol_deprotection_smarts)

        #Create mol for bromoacetic acid 
        self.bromoacetic_acid_smiles='OC(CBr)=O'
        self.bromoacetic_acid_mol=Chem.MolFromSmiles(self.bromoacetic_acid_smiles)


    def amine_displacement(self, amine:str, alkyl_bromide_intermediate:str):
        #create the mols
        self.alkyl_bromide_mol=Chem.MolFromSmiles(alkyl_bromide_intermediate)
        amine_mol=Chem.MolFromSmiles(amine)

        #Run the reaction
        product = self.displacement_reaction.RunReactants((amine_mol,self.alkyl_bromide_mol))[0][0]

        #Convert the product to smiles
        product_smiles = Chem.MolToSmiles(product)
        return(product_smiles)

    def amide_bond_formation1(self, amine_intermediate:str):
        #create the mol
        amine_mol = Chem.MolFromSmiles(amine_intermediate)

        #run the reaction
        product = self.amidation_reaction.RunReactants ([self.bromoacetic_acid_mol, amine_mol])

        #Convert the product to smiles
        product_smiles=Chem.MolToSmiles(product[0][0])
        return(product_smiles)


    def amide_bond_formation2(self, amine_intermediate:str):
        #create the mol
        amine_mol = Chem.MolFromSmiles(amine_intermediate)

        #block the existing amides from reacting as an amine
        amidep = Chem.MolFromSmarts('[N;$(NC=[O,S])]')
        for match in amine_mol.GetSubstructMatches(amidep):
            amine_mol.GetAtomWithIdx(match[0]).SetProp('_protected','1')
        
        #run the reaction
        product = self.amidation_reaction.RunReactants ([self.bromoacetic_acid_mol, amine_mol])

        #Convert the product to smiles
        product_smiles=Chem.MolToSmiles(product[0][0])
        return(product_smiles)

    def deprotect_peptoid(self, protected_peptoid_sequence:str, protected_peptoid_smiles:str):

        # Initialize the peptoid class object to for running deprotection reactions
        if self.int_deprotected_mol==None:
                self.int_deprotected_smiles=protected_peptoid_smiles
                self.int_deprotected_mol=Chem.MolFromSmiles(protected_peptoid_smiles)

        #Deprotect boc groups
        if any(match in protected_peptoid_sequence for match in ['G']):
            while self.gly_deprotection_reaction.RunReactantInPlace( self.int_deprotected_mol):
                pass
            Chem.SanitizeMol(self.int_deprotected_mol)
            self.int_deprotected_smiles=Chem.MolToSmiles(self.int_deprotected_mol)

        #Deprotect boc groups
        if any(match in protected_peptoid_sequence for match in ['D', 'K', 'R']):
            while self.boc_deprotection_reaction.RunReactantInPlace(self.int_deprotected_mol):
                pass
            Chem.SanitizeMol(self.int_deprotected_mol)
            self.int_deprotected_smiles=Chem.MolToSmiles(self.int_deprotected_mol)
        if any(match in protected_peptoid_sequence for match in ['W', 'H']): # The aromatic boc groups are handled separately 
            while self.Ar_boc_deprotection_reaction.RunReactantInPlace(self.int_deprotected_mol):
                pass
            Chem.SanitizeMol(self.int_deprotected_mol)
            self.int_deprotected_smiles=Chem.MolToSmiles(self.int_deprotected_mol)

        #Deprotect alcohols
        if any(match in protected_peptoid_sequence for match in ['S', 'Y']):
            while self.alcohol_deprotection_reaction.RunReactantInPlace(self.int_deprotected_mol):
                pass
            Chem.SanitizeMol(self.int_deprotected_mol)
            self.int_deprotected_smiles=Chem.MolToSmiles(self.int_deprotected_mol)
    
        #Deprotect carboxylates and obtain fully deprotected peptoid 
        while self.ester_deprotection_reaction.RunReactantInPlace(self.int_deprotected_mol):
            pass
        Chem.SanitizeMol(self.int_deprotected_mol)
        self.int_deprotected_smiles=Chem.MolToSmiles(self.int_deprotected_mol)

        return(self.int_deprotected_smiles)
