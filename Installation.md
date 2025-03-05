# Installation Guide

## Prerequisites
Before running the **SMILES to SDF Converter**, ensure you have:
✅ **Python 3.7+** installed.
✅ **RDKit** for molecular processing.
✅ **Pandas** for handling CSV files.

## Installation Steps
### 1️⃣ Install Python
Download and install Python from the [official Python website](https://www.python.org/downloads/).

### 2️⃣ Set Up a Virtual Environment (Optional but Recommended)
Create and activate a virtual environment:
```bash
# Create a virtual environment
python -m venv smiles_env

# Activate it (Windows)
smiles_env\Scripts\activate

# Activate it (Mac/Linux)
source smiles_env/bin/activate
```

### 3️⃣ Install Required Dependencies
Run the following command to install all necessary Python packages:
```bash
pip install pandas rdkit
```

If you're using **Anaconda**, install RDKit using:
```bash
conda install -c conda-forge rdkit
```

### 4️⃣ Verify the Installation
To check if RDKit is installed correctly, run:
```python
from rdkit import Chem
print(Chem.MolFromSmiles("CCO"))  # Should return a valid molecule object
```

### 5️⃣ Run the Script
Ensure your **CSV file** (`molecules.csv`) is ready and contains SMILES structures. Then, execute:
```bash
python main.py
```

### 6️⃣ Output File
After execution, check for the output file (`molecules_results.csv` or `.sdf`/`.json` based on your selection).

## Troubleshooting
### RDKit Installation Issues
If you get `No module named 'rdkit'`, use:
```bash
conda install -c conda-forge rdkit
```

### File Not Found Error
Make sure your **molecules.csv** file is in the correct directory.
```bash
ls  # Linux/Mac
dir # Windows
```

## Support
For any issues, contact **Sree Vasthav Upputoori** 

