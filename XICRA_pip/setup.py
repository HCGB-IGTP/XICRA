import setuptools
import glob

#######
def get_require_modules():
    """
    Get main python requirements modules
    """
    with open("./XICRA/config/python/python_requirement_summary.txt", 'r') as f:
        myModules = [line.strip().split(',')[0] for line in f]
    
    return myModules

#######
def get_version(file_VERSION):
    """
    Original code: PhiSpy setup.py 
    https://github.com/linsalrob/PhiSpy/blob/master/setup.py
    """
    with open(file_VERSION, 'r') as f:
        v = f.readline().strip()
    return v

long_description_text = ""
with open("README.md", "r") as fh:
    long_description_text = fh.read()

setuptools.setup(
    name="XICRA",
    version=get_version("XICRA/config/VERSION"),

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
    install_requires=get_require_modules(),

)
