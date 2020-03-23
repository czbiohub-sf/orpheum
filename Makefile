help:
	@echo "lint - check code style with flake8"
	@echo "test - run tests only"
	@echo "coverage - run tests and check code coverage"
	@echo "conda_install (recommended) - Install requirements "

test:
	py.test

coverage:
	coverage run --source khtools --omit="*/test*" --module py.test
	coverage report --show-missing

lint:
	# ignore:
	# E126 continuation line over-indented for hanging indent
	# E127 continuation line over-indented for visual indent
	flake8 --exclude docs,tests --ignore=E127,E126 khtools
	# stop the build if there are Python syntax errors or undefined names
    flake8 . --count --exclude docs,tests --select=E9,F63,F7,F82 --show-source --statistics
    # exit-zero treats all errors as warnings. The GitHub editor is 127 chars wide
    flake8 . --count --exclude docs,tests --exit-zero --max-complexity=10 --max-line-length=127 --statistics

conda_install:
	conda install --file conda_requirements.txt
	pip install -r requirements.txt
	pip install .

install:
	pip install -r requirements.txt
	pip install .
