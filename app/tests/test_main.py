from fastapi.testclient import TestClient


from main import app
from models import EXAMPLE


EXAMPLE_CHIRAL = {
   "reaction": "[C:1][C@H:2]([CH3:3])[I:4]>>[C:1][C@@H:2](Br)[CH3:3]",
   "reactants": [
       "CCCC[C@@H](I)C"
   ],
   "threshold": 10
}
EXAMPLE_MULTIPLE_REACTANTS = {
   "reaction": "[Br:1].[C:2]=[C:3]>>[Br:1][C:2][C:3]",
   "reactants": [
       "Br",
       "C=CC"
   ],
   "threshold": 10
}
EXAMPLE_NOT_MATCHING = {
   "reaction": "[C:1][C@H:2]([CH3:3])[I:4]>>[C:1][C@@H:2](Br)[CH3:3]",
   "reactants": [
       "CCCC[C@@H](F)C"
   ],
   "threshold": 10
}
EXAMPLE_INVALID_SMILES = {
   "reaction": "[C:1][C@H:2]([CH3:3])[I:4]>>[C:1][C@@H:2](Br)[CH3:3]",
   "reactants": [
       "blabla-CCCC[C@@H](F)C"
   ],
   "threshold": 10
}
EXAMPLE_INVALID_SMARTS = {
   "reaction": "blabla-[C:1][C@H:2]([CH3:3])[I:4]>>[C:1][C@@H:2](Br)[CH3:3]",
   "reactants": [
       "CCCC[C@@H](F)C"
   ],
   "threshold": 10
}




client = TestClient(app)




def test_home():
   response = client.get("/")
   assert response.status_code == 200
   assert "PerformReaction API" in response.text




def test_perform_reaction():
   response = client.post("/PerformReaction/", json=EXAMPLE)
   assert response.status_code == 200
   assert len(response.json()["rdkit_products"]) == 4
   assert len(response.json()["rdchiral_products"]) == 0




def test_perform_reaction_chiral():
   response = client.post('/PerformReaction/', json=EXAMPLE_CHIRAL)
   assert response.status_code == 200
   assert len(response.json()["rdkit_products"]) == 1
   assert len(response.json()["rdchiral_products"]) == 1




def test_perform_reaction_multiple_reactants():
   response = client.post('/PerformReaction/', json=EXAMPLE_MULTIPLE_REACTANTS)
   assert response.status_code == 200
   assert len(response.json()["rdkit_products"]) == 2
   assert len(response.json()["rdchiral_products"]) == 0




def test_perform_reaction_not_matching():
   response = client.post("/PerformReaction/", json=EXAMPLE_NOT_MATCHING)
   assert response.status_code == 200
   assert len(response.json()["rdkit_products"]) == 0
   assert len(response.json()["rdchiral_products"]) == 0




def test_perform_reaction_invalid_smiles():
   response = client.post("/PerformReaction/", json=EXAMPLE_INVALID_SMILES)
   assert response.status_code == 400
   assert response.json()["message"] == "Invalid SMILES"




def test_perform_reaction_invalid_smarts():
   response = client.post("/PerformReaction/", json=EXAMPLE_INVALID_SMARTS)
   assert response.status_code == 400
   assert response.json()["message"] == "Invalid SMARTS"
