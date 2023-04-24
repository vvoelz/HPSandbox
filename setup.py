
#!/usr/bin/env python
from setuptools import setup, find_packages
from os import path
from io import open
import sys
from pip._internal.req import parse_requirements

# parse_requirements() returns generator of pip.req.InstallRequirement objects
# We need two seperate requirement files due to the failure to install all at once.
requirements = parse_requirements('requirements.txt', session=False)
requirements = list(requirements)
try:
    requirements = [str(ir.req) for ir in requirements]
except:
    requirements = [str(ir.requirement) for ir in requirements]

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'README.md'), "r", encoding='utf-8') as f:
    long_description = f.read()


setup(
        name='hpsandbox',
        version='1.0.0b',
        #scripts=['hpsandbox'] ,
        author="Vincent Voelz",
        author_email="vvoelz@gmail.com",
        description="Python package to explore the 2D HP model of Chan and Dill.",
        long_description=long_description,
        long_description_content_type="text/markdown",
        url="https://github.com/vvoelz/HPSandbox",
        packages=find_packages(),
        install_requires=requirements,
        classifiers=[
            "Programming Language :: Python :: 3",
            "License :: OSI Approved :: MIT License",
            "Operating System :: OS Independent",
        ]
)





