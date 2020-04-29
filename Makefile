help:
	@echo "lint - check code style with flake8"
	@echo "test - run tests only"
	@echo "coverage - run tests and check code coverage"
	@echo "conda_install (recommended) - Install requirements "

test:
	py.test

coverage:
	coverage run --source sencha --omit="*/test*" --module py.test
	coverage report --show-missing

lint:
	# ignore:
	# E126 continuation line over-indented for hanging indent
	# E127 continuation line over-indented for visual indent
	flake8 --exclude docs,tests --ignore=E127,E126 sencha

conda_install:
	conda install --file conda_requirements.txt
	pip install -r requirements.txt
	pip install .

install:
	pip install -r requirements.txt
	pip install .

dist: FORCE
	python setup.py sdist

FORCE:
