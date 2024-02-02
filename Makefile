install:
	pip install --upgrade pip &&\
		python -m pip install -r requirements.txt

format:
	black *.py

lint:
	pylint --disable=R,C main.py

test: 
#pytest
	python -m pytest -vv --cov=main test/test_main.py 