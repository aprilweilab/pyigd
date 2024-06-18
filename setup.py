from setuptools import setup, find_packages

PACKAGE_NAME = "pyigd"
VERSION = "0.1"

setup(name=PACKAGE_NAME,
      packages=find_packages(),
      version=VERSION,
      description="Parser for IGD files",
      author="Drew DeHaas",
      author_email="",
      url="https://aprilweilab.github.io/",
      classifiers=[
          "Programming Language :: Python :: 3",
      ])
