#!/usr/bin/env python3
import setuptools
import glob
import HCGB.config.setup_module as setup_module

long_description_text = ""
with open("README.md", "r") as fh:
    long_description_text = fh.read()

setuptools.setup(
    name="XICRA",
    version=setup_module.get_version("./VERSION"),

    scripts=glob.glob('main/*'),
    author="Jose F. Sanchez-Herrero",

    author_email="jfbioinformatics@gmail.com",
    description="Small RNAseq pipeline for paired-end reads",

    long_description_content_type="text/markdown",
    long_description=long_description_text,
    url="https://github.com/HCGB-IGTP/XICRA/",
    packages=setuptools.find_packages(),
    license='MIT License',

    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    include_package_data=True,
    install_requires=setup_module.get_require_modules("XICRA/config/python/python_requirement_summary.txt"),

)
