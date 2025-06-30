import pandas as pd
import logging
import os
import glob
from rdkit import Chem
from rdkit.Chem import AllChem, Descriptors, rdmolfiles
from concurrent.futures import ThreadPoolExecutor, as_completed
import chardet  # For encoding detection

# Try to import Excel support libraries
try:
    import openpyxl
    EXCEL_SUPPORT = True
except ImportError:
    try:
        import xlrd
        EXCEL_SUPPORT = True
    except ImportError:
        EXCEL_SUPPORT = False
        print("  Warning: Excel support not available. Install openpyxl or xlrd for Excel files.")

# Configure logging
logging.basicConfig(
    filename="conversion.log", 
    level=logging.INFO, 
    format="%(asctime)s - %(levelname)s - %(message)s",
    filemode='w'  # Overwrite log file each run
)

def display_available_files():
    """Display available data files in the current directory"""
    # Support multiple file formats
    file_patterns = [
        "*.csv", "*.xlsx", "*.xls", "*.tsv", "*.txt", 
        "*.CSV", "*.XLSX", "*.XLS", "*.TSV", "*.TXT"
    ]
    
    all_files = []
    for pattern in file_patterns:
        all_files.extend(glob.glob(pattern))
    
    # Remove duplicates and sort
    all_files = sorted(list(set(all_files)))
    
    if not all_files:
        print(" No data files found in the current directory.")
        print("   Supported formats: CSV, Excel (.xlsx/.xls), TSV, TXT")
        return None
    
    print("\n Available files:")
    print("-" * 60)
    for i, file in enumerate(all_files, 1):
        try:
            file_size = os.path.getsize(file) / 1024  # Size in KB
            file_ext = os.path.splitext(file)[1].upper()
            print(f"{i:2d}. {file:<35} {file_ext:<6} ({file_size:.1f} KB)")
        except OSError:
            print(f"{i:2d}. {file:<35} ERROR  (Cannot read size)")
    
    return all_files

def select_file():
    """Allow user to select a file from available options"""
    files = display_available_files()
    if not files:
        return None
    
    while True:
        try:
            choice = input(f"\n Select file number (1-{len(files)}): ").strip()
            if not choice:
                print("Please enter a number.")
                continue
            
            file_index = int(choice) - 1
            if 0 <= file_index < len(files):
                selected_file = files[file_index]
                print(f" Selected: {selected_file}")
                return selected_file
            else:
                print(f" Please enter a number between 1 and {len(files)}")
        except ValueError:
            print(" Please enter a valid number.")
        except KeyboardInterrupt:
            print("\n Operation cancelled.")
            return None

def detect_encoding(file_path):
    """Detect the encoding of a text file"""
    import chardet
    
    try:
        with open(file_path, 'rb') as f:
            raw_data = f.read(10000)  # Read first 10KB for detection
        result = chardet.detect(raw_data)
        return result['encoding']
    except:
        return 'utf-8'  # Default fallback

def load_data_file(file_path):
    """Load data from CSV or Excel file with robust encoding handling"""
    try:
        file_ext = os.path.splitext(file_path)[1].lower()
        
        if file_ext == '.csv':
            # Try multiple encodings for CSV files
            encodings_to_try = ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1', 'utf-16']
            
            # First try to detect encoding
            detected_encoding = detect_encoding(file_path)
            if detected_encoding and detected_encoding not in encodings_to_try:
                encodings_to_try.insert(0, detected_encoding)
            
            df = None
            encoding_used = None
            
            for encoding in encodings_to_try:
                try:
                    print(f" Trying encoding: {encoding}")
                    df = pd.read_csv(file_path, encoding=encoding)
                    encoding_used = encoding
                    print(f" Successfully loaded with {encoding} encoding")
                    break
                except UnicodeDecodeError:
                    continue
                except Exception as e:
                    if "codec can't decode" in str(e):
                        continue
                    else:
                        raise e
            
            if df is None:
                # Last resort: try with errors='ignore'
                try:
                    print(" Trying with error handling...")
                    df = pd.read_csv(file_path, encoding='utf-8', errors='ignore')
                    encoding_used = 'utf-8 (with errors ignored)'
                    print(" Loaded with error handling")
                except Exception as e:
                    raise ValueError(f"Could not read CSV file with any encoding: {e}")
        
        elif file_ext in ['.xlsx', '.xls']:
            try:
                df = pd.read_excel(file_path, engine='openpyxl' if file_ext == '.xlsx' else 'xlrd')
            except Exception:
                # Try alternative engines
                try:
                    df = pd.read_excel(file_path)
                except Exception as e:
                    raise ValueError(f"Could not read Excel file: {e}")
        
        elif file_ext == '.tsv':
            # Handle tab-separated files
            encodings_to_try = ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1']
            df = None
            
            for encoding in encodings_to_try:
                try:
                    df = pd.read_csv(file_path, sep='\t', encoding=encoding)
                    print(f" TSV loaded with {encoding} encoding")
                    break
                except UnicodeDecodeError:
                    continue
        
        else:
            # Try to read as CSV anyway (sometimes files have wrong extensions)
            print(f"  Unknown extension {file_ext}, trying as CSV...")
            encodings_to_try = ['utf-8', 'latin-1', 'cp1252', 'iso-8859-1']
            df = None
            
            for encoding in encodings_to_try:
                try:
                    df = pd.read_csv(file_path, encoding=encoding)
                    print(f" File loaded as CSV with {encoding} encoding")
                    break
                except:
                    continue
            
            if df is None:
                raise ValueError(f"Unsupported file format: {file_ext}")
        
        if df is None or df.empty:
            raise ValueError("The input file is empty or could not be read.")
        
        # Clean up column names (remove extra spaces, special characters)
        df.columns = df.columns.str.strip()
        
        print(f" Loaded {len(df)} rows from {file_path}")
        if encoding_used:
            print(f" Encoding used: {encoding_used}")
        print(f" Columns found: {', '.join(df.columns.tolist())}")
        
        # Show first few rows for verification
        print(f"\n Preview of data:")
        print(df.head(3).to_string())
        
        return df
    
    except Exception as e:
        logging.error(f"Error reading file {file_path}: {e}")
        print(f" Error reading file: {e}")
        return None

def identify_columns(df):
    """Identify Name and SMILES columns with user confirmation"""
    columns = df.columns.tolist()
    
    # Try to auto-detect common column names
    name_col = None
    smiles_col = None
    
    # Common name column variations
    name_variants = ['name', 'compound', 'molecule', 'title', 'id', 'compound_name', 'mol_name']
    for col in columns:
        if col.lower() in name_variants or 'name' in col.lower():
            name_col = col
            break
    
    # Common SMILES column variations
    smiles_variants = ['smiles', 'smile', 'canonical_smiles', 'smiles_string']
    for col in columns:
        if col.lower() in smiles_variants or 'smiles' in col.lower():
            smiles_col = col
            break
    
    print(f"\n Auto-detected columns:")
    print(f"   Name column: {name_col if name_col else 'Not found'}")
    print(f"   SMILES column: {smiles_col if smiles_col else 'Not found'}")
    
    # Let user confirm or select manually
    print(f"\n Available columns:")
    for i, col in enumerate(columns, 1):
        sample_data = str(df[col].dropna().iloc[0])[:50] if not df[col].dropna().empty else "No data"
        print(f"{i:2d}. {col:<20} (Sample: {sample_data})")
    
    # Confirm name column
    if name_col:
        confirm = input(f"\n Use '{name_col}' as Name column? (y/n): ").strip().lower()
        if confirm != 'y':
            name_col = None
    
    if not name_col:
        while True:
            try:
                choice = input(f"  Select Name column number (1-{len(columns)}): ").strip()
                idx = int(choice) - 1
                if 0 <= idx < len(columns):
                    name_col = columns[idx]
                    break
                else:
                    print(f"Please enter a number between 1 and {len(columns)}")
            except ValueError:
                print("Please enter a valid number.")
    
    # Confirm SMILES column
    if smiles_col:
        confirm = input(f"\n Use '{smiles_col}' as SMILES column? (y/n): ").strip().lower()
        if confirm != 'y':
            smiles_col = None
    
    if not smiles_col:
        while True:
            try:
                choice = input(f"🧬 Select SMILES column number (1-{len(columns)}): ").strip()
                idx = int(choice) - 1
                if 0 <= idx < len(columns):
                    smiles_col = columns[idx]
                    break
                else:
                    print(f"Please enter a number between 1 and {len(columns)}")
            except ValueError:
                print("Please enter a valid number.")
    
    return name_col, smiles_col

# Global variable for hydrogen handling mode
HYDROGEN_MODE = 3  # Default: Smart mode

def process_molecule_with_mode(data):
    """Process molecule based on selected hydrogen handling mode"""
    idx, name, smiles = data
    
    # Handle missing or invalid data
    if pd.isna(name) or not name:
        name = f"Molecule_{idx}"
    
    if pd.isna(smiles) or not smiles:
        logging.warning(f"Entry {idx}: Missing SMILES for {name}")
        return None
    
    # Clean SMILES string
    smiles = str(smiles).strip()
    
    try:
        # Create molecule from SMILES
        mol = Chem.MolFromSmiles(smiles)
        if mol is None:
            logging.warning(f"Entry {idx}: Invalid SMILES '{smiles}' for {name}")
            return None
        
        if HYDROGEN_MODE == 1:
            # Mode 1: Preserve original structure (no hydrogen addition)
            final_mol = mol
            final_mol.SetProp("Structure", "Original")
            final_mol.SetProp("Hydrogens", "Implicit")
            logging.info(f"Entry {idx}: Preserved original structure for {name}")
            
        elif HYDROGEN_MODE == 2:
            # Mode 2: Always add explicit hydrogens
            mol_with_h = Chem.AddHs(mol)
            embed_result = AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG())
            
            if embed_result == 0:
                try:
                    AllChem.MMFFOptimizeMolecule(mol_with_h, maxIters=200)
                except:
                    pass
                final_mol = mol_with_h
                final_mol.SetProp("Structure", "3D_with_H")
            else:
                final_mol = mol_with_h
                final_mol.SetProp("Structure", "2D_with_H")
            
            final_mol.SetProp("Hydrogens", "Explicit")
            logging.info(f"Entry {idx}: Added explicit hydrogens for {name}")
            
        else:
            # Mode 3: Smart mode (default)
            mol_with_h = Chem.AddHs(mol)
            embed_result = AllChem.EmbedMolecule(mol_with_h, AllChem.ETKDG())
            
            if embed_result == 0:
                # 3D successful - keep hydrogens and optimize
                try:
                    AllChem.MMFFOptimizeMolecule(mol_with_h, maxIters=200)
                except:
                    pass
                final_mol = mol_with_h
                final_mol.SetProp("Structure", "3D")
                final_mol.SetProp("Hydrogens", "Added_for_3D")
                logging.info(f"Entry {idx}: Generated 3D coordinates with H for {name}")
            else:
                # 3D failed - use original structure without added hydrogens
                final_mol = mol
                final_mol.SetProp("Structure", "2D_Original")
                final_mol.SetProp("Hydrogens", "Implicit")
                logging.info(f"Entry {idx}: Kept original 2D structure for {name}")
        
        # Set molecule properties
        final_mol.SetProp("_Name", str(name))
        final_mol.SetProp("ID", str(idx))
        final_mol.SetProp("SMILES", smiles)
        
        # Calculate properties using original molecule (more accurate)
        mol_weight = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Descriptors.NumHDonors(mol)
        hba = Descriptors.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        rotatable_bonds = Descriptors.NumRotatableBonds(mol)
        
        # Set molecular properties
        final_mol.SetProp("MolWt", f"{mol_weight:.2f}")
        final_mol.SetProp("LogP", f"{logp:.2f}")
        final_mol.SetProp("HBD", str(hbd))
        final_mol.SetProp("HBA", str(hba))
        final_mol.SetProp("TPSA", f"{tpsa:.2f}")
        final_mol.SetProp("RotBonds", str(rotatable_bonds))
        
        return final_mol
        
    except Exception as e:
        logging.error(f"Entry {idx}: Error processing {name} with SMILES '{smiles}': {e}")
        return None

def convert_to_sdf(input_file):
    """Convert CSV/Excel file to SDF format with structure preservation"""
    print("\n Starting conversion process...")
    
    # Load the data file
    df = load_data_file(input_file)
    if df is None:
        return
    
    # Identify the name and SMILES columns
    name_col, smiles_col = identify_columns(df)
    
    # Ask user about hydrogen handling preference
    print(f"\n Hydrogen Handling Options:")
    print("1. Preserve original structure (recommended for visualization)")
    print("2. Add explicit hydrogens for 3D coordinates")
    print("3. Smart mode (add H only if 3D generation succeeds)")
    
    while True:
        try:
            h_choice = input("Select option (1-3) [default: 3]: ").strip()
            if not h_choice:
                h_choice = "3"
            if h_choice in ['1', '2', '3']:
                break
            print("Please enter 1, 2, or 3")
        except KeyboardInterrupt:
            print("\nOperation cancelled.")
            return
    
    global HYDROGEN_MODE
    HYDROGEN_MODE = int(h_choice)
    
    # Prepare data for processing
    df_clean = df[[name_col, smiles_col]].copy()
    df_clean = df_clean.dropna(subset=[smiles_col])  # Remove rows with missing SMILES
    
    print(f"\n Processing {len(df_clean)} molecules...")
    print(f" Hydrogen mode: {['', 'Preserve original', 'Add explicit H', 'Smart mode'][HYDROGEN_MODE]}")
    
    # Create output filename
    base_name = os.path.splitext(input_file)[0]
    output_sdf = f"{base_name}_converted.sdf"
    
    # Prepare data for parallel processing
    molecules_data = [(i+1, row[name_col], row[smiles_col]) 
                     for i, (_, row) in enumerate(df_clean.iterrows())]
    
    # Process molecules in parallel
    successful_molecules = []
    failed_count = 0
    
    print("  Processing molecules...")
    with ThreadPoolExecutor(max_workers=4) as executor:
        # Submit all tasks
        future_to_data = {executor.submit(process_molecule_with_mode, data): data for data in molecules_data}
        
        # Collect results as they complete
        for future in as_completed(future_to_data):
            mol = future.result()
            if mol is not None:
                successful_molecules.append(mol)
            else:
                failed_count += 1
            
            # Show progress
            completed = len(successful_molecules) + failed_count
            if completed % 10 == 0 or completed == len(molecules_data):
                print(f"   Progress: {completed}/{len(molecules_data)} molecules processed")
    
    # Write to SDF file
    if successful_molecules:
        try:
            with Chem.SDWriter(output_sdf) as writer:
                for mol in successful_molecules:
                    writer.write(mol)
            
            print(f"\n Successfully converted {len(successful_molecules)} molecules")
            print(f" Output saved as: {output_sdf}")
            
            if failed_count > 0:
                print(f"  {failed_count} molecules failed to convert (check conversion.log for details)")
            
            # Display summary with structure information
            structure_3d = sum(1 for mol in successful_molecules if mol.HasProp("Structure") and mol.GetProp("Structure") == "3D")
            structure_2d = len(successful_molecules) - structure_3d
            
            print(f"\n Conversion Summary:")
            print(f"   Input file: {input_file}")
            print(f"   Total molecules: {len(molecules_data)}")
            print(f"   Successfully converted: {len(successful_molecules)}")
            print(f"   With 3D coordinates: {structure_3d}")
            print(f"   With 2D structure: {structure_2d}")
            print(f"   Failed conversions: {failed_count}")
            print(f"   Output file: {output_sdf}")
            
        except Exception as e:
            logging.error(f"Error writing SDF file: {e}")
            print(f" Error writing SDF file: {e}")
    else:
        print(" No molecules were successfully converted!")

def main():
    """Main function to run the converter"""
    print(" Molecule Converter - Multiple Formats to SDF")
    print("=" * 60)
    print(" Supported formats: CSV, Excel (.xlsx/.xls), TSV, TXT")
    print(" Auto-detects file encoding (UTF-8, Latin-1, CP1252, etc.)")
    
    if not EXCEL_SUPPORT:
        print("  Note: Install 'openpyxl' for better Excel support")
    
    try:
        # Select input file
        input_file = select_file()
        if not input_file:
            return
        
        # Convert to SDF
        convert_to_sdf(input_file)
        
    except KeyboardInterrupt:
        print("\n Conversion cancelled by user.")
    except Exception as e:
        logging.error(f"Unexpected error in main: {e}")
        print(f" An unexpected error occurred: {e}")
        print(" Try installing missing dependencies: pip install chardet openpyxl xlrd")

if __name__ == "__main__":
    main()
