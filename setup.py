# coding: utf-8

from setuptools import find_packages, setup
from pathlib import Path

# The directory containing this file
HERE = Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

setup(
    name='pseudogene',
    version='1.0.0.dev3',
    packages=find_packages(),
    url='https://github.com/felix-yczhong/pseudogene.git',
    license='GNU General Public License v.3.0',
    author='Yu-Chang Zhong',
    author_email='r09922111@ntu.edu.tw',
    long_description=README,
    long_description_content_type="text/markdown",
    python_requires='>=3.8',
    install_requires=['click', 'cython', 'numpy', 'pysam', 'joblib', 'openpyxl', 'pandas', 'smaca==1.2.3', 'gffutils'],
    include_package_data=True,
    entry_points={'console_scripts': [
                    'pseudogene = pseudogene.main:main']})