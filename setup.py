# coding: utf-8

from setuptools import find_packages, setup
from pathlib import Path

# The directory containing this file
HERE = Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
    name='pseudogene',
    version='1.0.0.dev1',
    packages=find_packages(),
    url='',
    license='GNU General Public License v.3.0',
    author='Yu-Chang Zhong',
    author_email=
    '',
    long_description=README,
    long_description_content_type="text/markdown",
    python_requires='>=3.8',
    install_requires=['click', 'cython', 'numpy', 'pysam', 'joblib', 'openpyxl', 'pandas', 'smaca==1.2.3'],
    include_package_data=True,
    entry_points={'console_scripts': ['pseudogene = pseudogene.main:main']})