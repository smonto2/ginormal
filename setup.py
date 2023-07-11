import setuptools
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setuptools.setup(
    name='gin',
    version='0.0.1',
    packages=setuptools.find_packages(),
    author='Santiago Montoya-Blandón, Zhilang Xia, Cheng Ding, Juan Estrada',
    author_email='Santiago.Montoya-Blandon@glasgow.ac.uk',
    description='Generalized Inverse Normal (GIN) distribution density and generation',
    long_description=long_description,           #in README file
    long_description_content_type='text/markdown',
    install_requires=['numpy','scipy','cardano_method'],           #the basis package needed
    url='https://github.com/smonto2/GIN',
    classifiers=[
        "Programming Language :: Python :: 3.8", #"Programming Language :: Python :: 3", "Programming Language" represents the top-level category, "Python" represents the second-level category, and "3.8" represents the third-level category.
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",    #the software package is not dependent on any specific operating system
    ],
    python_requires='>=3.8',
)
