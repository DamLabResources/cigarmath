VENV_NAME?=venv
VENV_ACTIVATE=. $(VENV_NAME)/bin/activate
PYTHON=${VENV_NAME}/bin/python3

.PHONY: venv

venv: $(VENV_NAME)/bin/activate

$(VENV_NAME)/bin/activate: requirements.txt
	test -d $(VENV_NAME) || python3 -m venv $(VENV_NAME)
	${PYTHON} -m pip install --upgrade pip
	${PYTHON} -m pip install --upgrade setuptools wheel
	. $(VENV_NAME)/bin/activate; pip install -r requirements.txt
	${PYTHON} setup.py install
	touch $(VENV_NAME)/bin/activate

update_env:
	${PYTHON} setup.py install
	touch $(VENV_NAME)/bin/activate
    
clean:
	rm -rf $(VENV_NAME)

test: venv
	${PYTHON} -m pytest

run: venv
	${PYTHON} your_script.py

nbstrip:
	nbstripout --keep-output notebooks/*.ipynb
