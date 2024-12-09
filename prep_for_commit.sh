#!/bin/bash
set -ev

flake8 pyigd --count --select=E9,F63,F7,F82,F401 --show-source --statistics

black pyigd/ setup.py test/ examples/ --check 

mypy pyigd --no-namespace-packages --ignore-missing-imports

pytest test/

