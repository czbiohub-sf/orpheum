PYTHON ?= python

help:
	@echo "lint - check code style with flake8"
	@echo "code_format - check code style with black"
	@echo "test - run tests only"
	@echo "coverage - run tests and check code coverage"
	@echo "conda_install (recommended) - Install requirements "

clean:
	$(PYTHON) setup.py clean --all

install:
	pip install -r requirements.txt
	pip install .

dist: FORCE
	$(PYTHON) setup.py sdist

test:
	pip install -e '.[test]'
	$(PYTHON) -m pytest

coverage:
	$(PYTHON) setup.py clean --all
	$(PYTHON) -m pytest --cov=. --cov-report term-missing

FORCE:

code_format:
	black .

lint:
	# ignore:
	# E126 continuation line over-indented for hanging indent
	# E127 continuation line over-indented for visual indent
	flake8 --exclude docs,tests --ignore=E127,E126 orpheum

conda_install:
	conda install --file conda_requirements.txt
	pip install -r requirements.txt
	pip install .

