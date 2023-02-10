import logging
import rdkit.Chem as Chem
from fastapi import FastAPI, Request
from fastapi.responses import HTMLResponse, JSONResponse
from rdchiral.main import rdchiralRunText


from models import ReactionRequest, ReactionResponse
from validators import validate_reactants, validate_reaction


app = FastAPI(title="Run SMILES through reaction SMARTS")




@app.get("/", response_class=HTMLResponse)
def home():
   return '<h1>PerformReaction API</h1>\n<a href="docs">[Documentation link]</a>'




@app.post("/PerformReaction/")
def perform_reaction(data: ReactionRequest) -> ReactionResponse:
   mols = validate_reactants(data.reactants)
   rxn = validate_reaction(data.reaction)


   # Run reactants (rdkit)
   rdkit_products = {tuple(sorted([Chem.MolToSmiles(mol) for mol in mols]))
                     for mols in rxn.RunReactants(mols, maxProducts=data.threshold)}


   # Run reactants (rdchiral)
   try:
       if len(data.reactants) == 1:
           rdchiral_products = rdchiralRunText(data.reaction, data.reactants[0])[0].split('.')
       else:
           logging.error(f'rdchiral error: more than 1 reactant')
           rdchiral_products = []
   except Exception as e:
       logging.error(f'rdchiral error: {e}')
       rdchiral_products = []


   return ReactionResponse(reaction=data.reaction, reactants=data.reactants,
                           rdkit_products=rdkit_products, rdchiral_products=rdchiral_products)




@app.exception_handler(ValueError)
async def value_error_handler(request: Request, exc: ValueError):
   return JSONResponse(
       status_code=400,
       content={"message": str(exc)},
   )
