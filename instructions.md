Let's begin by installing any and all requirements using docker


```bash
cd app
docker build . -f Dockerfile.txt
```
Paste this into your terminal to load up the 
```bash
uvicorn main:app --reload 
```

Run the test sent in the e-mail
```bash
curl -X POST http://127.0.0.1:8000/PerformReaction/ \
  -H "Content-Type: application/json" \
  -d @test.json
```
Run a test to show unique products in a \
reaction with alcohols in a list
```bash
curl -X POST http://127.0.0.1:8000/PerformReaction/ \
  -H "Content-Type: application/json" \
  -d @test_alcohols.json
```
Test to see how invalid SMILES work
```bash
curl -X POST http://127.0.0.1:8000/PerformReaction/ \
  -H "Content-Type: application/json" \
  -d @test_badSMILES.json
```
Test to see how invalid SMARTS work
```bash
curl -X POST http://127.0.0.1:8000/PerformReaction/ \
  -H "Content-Type: application/json" \
  -d @test_badSMARTS.json
```
