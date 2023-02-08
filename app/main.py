from fastapi import FastAPI, Request
from pydantic import BaseModel
from fastapi.responses import JSONResponse
import rdkit.Chem as Chem
import rdkit.Chem.AllChem as AllChem
import rdchiral
from rdchiral.main import rdchiralRunText
from typing import List
from itertools import chain
import pandas as pd


app = FastAPI(title="Run SMILES through reaction SMARTS")


class ReactionHandler(BaseModel):
    reaction: str
    reactants: List[str]
    threshold: int


@app.get("/")
def home():
    """Documentation link"""
    return "The documentation can be found at http://127.0.0.1:8000/docs"


@app.exception_handler(ValueError)
async def errors(request: Request, exc: ValueError):
    return JSONResponse(
        status_code=400,
        content={"message": str(exc)},
    )


@app.post("/PerformReaction/")
def perform_reaction(data: ReactionHandler):
    reactants = data.reactants
    reaction = data.reaction
    thresh = data.threshold
    mols = [Chem.MolFromSmiles(smile) for smile in reactants]

    for molecule in mols:
        if molecule is None:
            raise ValueError('Invalid SMILES')

    rxn = AllChem.ReactionFromSmarts(reaction)
    if not rxn:
        raise ValueError('Invalid SMARTS')

    products = pd.DataFrame()
    products['reactants'] = reactants
    products['rdkit_products'] = [rxn.RunReactants((mol, ), maxProducts=thresh) for mol in mols]

    all_products = []
    for index, row in products.iterrows():
        all_products.append(list(chain.from_iterable(products['rdkit_products'].loc[index])))

    smiles = []
    for index, row in enumerate(all_products):
        molecules = all_products[index]
        smi = [Chem.MolToSmiles(mol) for mol in molecules]
        if smi not in smiles:
            smiles.append(smi)

    unique = []
    result = []
    for i in smiles:
        temp = []
        for x in i:
            if x not in unique:
                unique.append(x)
                temp.append(x)
        result.append(temp)

    products['rdkit_productSMILES'] = result

    try:
        products['rdchiral_productSMILES'] = [rdchiral.main.rdchiralRunText(reaction, mol) for mol in reactants]
    except Exception as e:
        print("Error occurred while performing rdchiralRunText: {}".format(str(e)))

    products = products.drop(columns='rdkit_products')

    return products.to_json()
