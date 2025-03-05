# SMILES to SDF Converter

## Overview
This project provides a Python script that converts **SMILES (Simplified Molecular Input Line Entry System)** strings from a CSV file into **SDF (Structure Data File)** format, while also calculating molecular properties like **Molecular Weight (MolWt), LogP, Hydrogen Bond Donors (HBD), and Hydrogen Bond Acceptors (HBA)** using RDKit.

## Features
✅ Converts **SMILES to SDF** format.
✅ Supports **CSV, JSON, and SDF** output formats.
✅ Uses **Parallel Processing** for speed improvement.
✅ Includes **Molecular Properties Calculation** (MolWt, LogP, HBD, HBA).
✅ **Error Handling & Logging** for invalid SMILES or missing data.

## Installation
### Prerequisites
Ensure you have **Python 3.7+** installed.

### Install Dependencies
```bash
pip install pandas rdkit
```

If you're using **Anaconda**, install RDKit via:
```bash
conda install -c conda-forge rdkit
```

## Usage
### Running the Script
Ensure you have a CSV file (e.g., `molecules.csv`) with the following format:
```
Name,SMILES
Aspirin,CC(=O)Oc1ccccc1C(=O)O
Paracetamol,CC(=O)Nc1ccc(O)cc1
...
```
Then, run:
```bash
python main.py
```

### Selecting Output Format
The script supports CSV, JSON, and SDF outputs. To specify the format:
```python
convert_smiles_csv("molecules.csv", "output.csv", output_format="sdf")
```
Possible formats:
- `csv` → Stores results in CSV format
- `json` → Saves output in JSON format
- `sdf` → Generates an SDF file

## Troubleshooting
### RDKit Import Error
If you see `No module named 'rdkit'`, install RDKit using:
```bash
conda install -c conda-forge rdkit
```

### File Not Found Error
Ensure the **CSV file exists** in the same directory as `main.py`. You can check using:
```bash
ls  # Linux/Mac
dir # Windows
```

## License
This project is open-source under the **MIT License**.

## Author
Developed by **Sree Vasthav Upputoori** - https://github.com/SreeVasthav-Upputoori 
