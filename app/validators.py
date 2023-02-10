from typing import List
from rdkit import Chem
from rdkit.Chem import AllChem

def validate_reactants(reactants: List[str]) -> List[Chem.Mol]:
   mols = []
   for reactant in reactants:
       mol = Chem.MolFromSmiles(reactant)
       if mol is None:
           raise ValueError('Invalid SMILES')
       mols.append(mol)
   return mols



def validate_reaction(reaction: str) -> AllChem.ChemicalReaction:
   try:
       rxn = AllChem.ReactionFromSmarts(reaction)
   except ValueError:
       raise ValueError('Invalid SMARTS')
   return rxn
