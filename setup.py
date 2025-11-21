from setuptools import setup, find_packages
import os

PACKAGE_NAME = "pyigd"
VERSION = "1.4"

THISDIR = os.path.dirname(os.path.realpath(__file__))
with open(os.path.join(THISDIR, "README.md")) as f:
    long_description = f.read()

setup(
    name=PACKAGE_NAME,
    packages=find_packages(),
    version=VERSION,
    description="Parser and writer of IGD files for genetic data",
    author="Drew DeHaas",
    author_email="",
    url="https://aprilweilab.github.io/",
    classifiers=[
        "Programming Language :: Python :: 3",
    ],
    long_description=long_description,
    long_description_content_type="text/markdown",
    license="MIT",
    keywords=[
        "population genetics",
        "popgen",
        "statistical genetics",
        "statgen",
        "computational biology",
        "compbio",
    ],
)
