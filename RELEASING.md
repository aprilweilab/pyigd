# Generating a new release

1. Bump version number in `setup.py`
2. Build packages and push to PyPi:
3. Update `conda/meta.yaml` with the latest version number and sha256 hash
4. Commit, push, and merge to main
5. Tag main with the version `vX.Y` and push

## Building the packages

```
rm -rf dist
python setup.py sdist
python setup.py bdist_wheel
python3 -m twine upload dist/*
```
