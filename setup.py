from setuptools import setup, find_packages

setup(
    name="csv_to_sdf_converter",
    version="1.0.0",
    author="Sree Vasthav Upputoori",
    author_email="sreevasthav.upputoori@gmail.com",
    description="A tool to convert SMILES structures to SDF format and compute molecular properties using RDKit.",
    long_description=open("README.md").read(),
    long_description_content_type="text/markdown",
    url="https://github.com/your-repo/smiles-to-sdf",  # Update with actual GitHub repo
    packages=find_packages(),
    install_requires=[
        "pandas",
        "rdkit-pypi"
    ],
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires='>=3.7',
    entry_points={
        "console_scripts": [
            "smiles2sdf=main:convert_smiles_csv"
        ]
    },
)
