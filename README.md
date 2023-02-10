# PerformReaction
Microservice running reactants with rdkit and rdchiral.

## Installation, testing and serving
### With docker (recommended for users)
```bash
docker build -t perform_reaction .
docker run -p 80:80 perform_reaction
```
### With venv (recommended for developers)
```bash
python3.11 -m venv env
source env/bin/activate
pip install -r requirements.txt
export PYTHONPATH=$PYTHONPATH:$PWD/app
cd app
pytest tests
uvicorn main:app --host 0.0.0.0 --port 80 --reload
```