from typing import List, Set, Tuple
from pydantic import BaseModel


EXAMPLE = {
   "reaction": "[c:8]-[c:6]>>[c:8][I:55].[B:99][c:6]",
   "reactants": [
       "CC1=CC=C(C=C1)C1=CC(=CC=C1C)C1=CC(C)=CC(C)=C1"
   ],
   "threshold": 10
}

class ReactionBase(BaseModel):
   reaction: str
   reactants: List[str]

class ReactionRequest(ReactionBase):
   threshold: int

   class Config:
       schema_extra = {"example": EXAMPLE}

class ReactionResponse(ReactionBase):
   rdkit_products: Set[Tuple[str, ...]]
   rdchiral_products: List[str]
