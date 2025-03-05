import pandas as pd
import logging
import json
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdmolfiles
from concurrent.futures import ThreadPoolExecutor

# Configure logging
logging.basicConfig(filename="conversion.log", level=logging.INFO, format="%(asctime)s - %(levelname)s - %(message)s")

def process_molecule(data):
    if len(data) != 3:  # Ensure correct tuple size
        logging.warning(f"Skipping malformed entry: {data}")
        return [None, None, None, "Invalid Entry", None, None, None, None]
    
    idx, name, smiles = data
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        logging.warning(f"Skipping invalid SMILES: {smiles}")
        return [idx, name, smiles, "Invalid SMILES", None, None, None, None]
    
    mol = Chem.AddHs(mol)
    AllChem.EmbedMolecule(mol, AllChem.ETKDG())
    sdf_str = rdmolfiles.MolToMolBlock(mol)
    
    # Calculate molecular properties
    mol_weight = Descriptors.MolWt(mol)
    logp = Descriptors.MolLogP(mol)
    hbd = Descriptors.NumHDonors(mol)
    hba = Descriptors.NumHAcceptors(mol)
    
    return [idx, name, smiles, sdf_str, mol_weight, logp, hbd, hba]

def convert_smiles_csv(input_csv, output_file, output_format="csv"):
    try:
        df = pd.read_csv(input_csv)
        if df.empty:
            logging.error("The input CSV file is empty.")
            return
    except Exception as e:
        logging.error(f"Error reading CSV: {e}")
        return
    
    if "Name" not in df.columns or "SMILES" not in df.columns:
        logging.error("The CSV file must contain 'Name' and 'SMILES' columns.")
        return
    
    smiles_list = df[["Name", "SMILES"]].dropna().values.tolist()
    print("First 5 entries:", smiles_list[:5])  # Debugging print
    
    with ThreadPoolExecutor(max_workers=4) as executor:
        results = list(executor.map(process_molecule, [(idx, name, smiles) for idx, (name, smiles) in enumerate(smiles_list, start=1)]))
    
    if output_format == "csv":
        result_df = pd.DataFrame(results, columns=["S.No.", "Name", "SMILES", "SDF", "MolWt", "LogP", "HBD", "HBA"])
        result_df.to_csv(output_file, index=False)
    elif output_format == "json":
        with open(output_file.replace(".csv", ".json"), "w") as f:
            json.dump(results, f, indent=4)
    elif output_format == "sdf":
        with Chem.SDWriter(output_file.replace(".csv", ".sdf")) as writer:
            for row in results:
                mol = Chem.MolFromSmiles(row[2])
                if mol:
                    writer.write(mol)
    
    logging.info(f"✅ Results saved in: {output_file}")
    print(f"✅ Results saved in: {output_file}")

# Example Usage
if __name__ == "__main__":
    input_csv_file = "molecules.csv"
    output_csv_file = "molecules_results.csv"
    convert_smiles_csv(input_csv_file, output_csv_file, output_format="csv")
