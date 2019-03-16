"""A setuptools based setup module.
See:
https://packaging.python.org/en/latest/distributing.html
https://github.com/pypa/sampleproject
"""

# Always prefer setuptools over distutils
from setuptools import setup

# To use a consistent encoding
from codecs import open
from os import path

here = path.abspath(path.dirname(__file__))

# Get the long description from the README file
with open(path.join(here, 'readme.md'), encoding='utf-8') as f:
    long_description = f.read()

# Arguments marked as "Required" below must be included for upload to PyPI.
# Fields marked as "Optional" may be commented out.

setup(
    name='read_lif',
    version='0.1.0',
    description='A Python module for loading lif file as numpy array',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/yangyushi/read_lif',
    author='Yushi Yang',
    author_email='yangyushi1992@icloud.com',  # Optional
    py_modules=["read_lif"],
    install_requires=['numpy'],
    python_requires='>=2.5'
)
