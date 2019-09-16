
PYTHON=./env/bin/python
CONDA=conda

all: env

setup: env pip

env:
	${CONDA} env create -f requirements.yml -p env

pip: env
	${PYTHON} -m pip install -r requirements.txt --no-cache-dir


test:
	${PYTHON} -m pytest test.py
