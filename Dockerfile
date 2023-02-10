FROM python:3.11

COPY requirements.txt /root
COPY app /root/app

RUN pip install -r /root/requirements.txt

ENV PYTHONPATH "$PYTHONPATH:/root/app"
WORKDIR /root/app

RUN pytest tests

CMD ["uvicorn", "main:app", "--host", "0.0.0.0", "--port", "80"]
