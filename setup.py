from setuptools import setup
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
    version='0.2.1',
    description='A Python module for loading lif file as numpy array',
    long_description=long_description,
    long_description_content_type='text/markdown',
    url='https://github.com/yangyushi/read_lif',
    author='Yushi Yang',
    author_email='yangyushi1992@icloud.com',  # Optional
    packages=["read_lif"],
    package_dir={'read_lif': 'read_lif'},
    install_requires=['numpy'],
    python_requires='>=2.5'
)
