import setuptools
import glob

long_description_text = ""
with open("README.md", "r") as fh:
    long_description_text = fh.read()

setuptools.setup(
    name="XICRA",
    version="0.9.2.2",

    scripts=glob.glob('main/*'),
    author="Jose F. Sanchez-Herrero",

    author_email="jfbioinformatics@gmail.com",
    description="",

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

    install_requires=[
        'pandas', 'patool', 'termcolor', 'cutadapt', 'mirtop', 'pysam', 'pybedtools', 'biopython', 'multiqc', 'HCGB'
    ],
)
