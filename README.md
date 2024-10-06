# PeptoidGeneration
Generate peptoid SMILES from amine building blocks (monomers).

See the incuded noteboof for an example of intended use.
* The amine monomers are input as a dictionary of smiles keyed by single-letter code for each monomer.
* The generator accepts a string representing the sequence of the desired peptoid.
* Peptoid generation is performed using iterative synthesis cycles akin to what is done in real-life:
    1. Amine Displacement: The amine monomer displaces a haloacetic acid to yield an a-amino substituted intermediate.
    2. Amide Coupling: The a-amino intermediate is coupled with the next unit of haloacetic acid.

To prevent cross reactivity, the smiles of the monomers need to protected:
* Amines not the intended as the reaction center need to be Boc protected.
* Carboxylic acids must be protected as the t-butyl ester. 

The generator will perform a final global deprotection of these groups.
