#!/usr/bin/env bash
python3 -m venv XICRA_env

source XICRA_env

pip install --upgrade pip

pip install configparser

pip install numpy
pip install cython

pip install biopython

pip install pandas
pip install xtract
