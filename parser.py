import pandas as pd

# Essential Fields Required for Quality Control
# These columns must be present in the dataset.
# They are also used to remove incomplete rows.

ESSENTIAL_FIELDS = [
    "Plasmid_Variant_Index",
    "Parent_Plasmid_Variant",
    "Directed_Evolution_Generation",
    "Assembled_DNA_Sequence",
    "DNA_Quantification_fg",
    "Protein_Quantification_pg"
]

# File Parsing Function

def parse_file(file_path):
    """
    Load TSV or JSON file and return pandas DataFrame.
    """

    # Check file extension and load accordingly
    if file_path.endswith(".tsv"):
        df = pd.read_csv(file_path, sep="\t")

    elif file_path.endswith(".json"):
        df = pd.read_json(file_path)

    else:
        # Raise error if format is not supported
        raise ValueError("Unsupported file format")

    return df

# Quality Control Function

def run_quality_control(df):
    """
    Validate schema and remove incomplete rows.
    
    Steps:
    1. Check if essential columns exist.
    2. Remove rows with missing values in essential fields.
    3. Separate clean and rejected rows.
    """

    # Check for missing essential columns
    missing_columns = [col for col in ESSENTIAL_FIELDS if col not in df.columns]

    if missing_columns:
        raise ValueError(f"Missing essential columns: {missing_columns}")

    # Remove rows with missing values in essential fields
    clean_df = df.dropna(subset=ESSENTIAL_FIELDS)

    # Identify rejected rows (rows removed during QC)
    rejected_df = df.loc[~df.index.isin(clean_df.index)]

    return clean_df, rejected_df

# Combined Parsing + QC Function

def parse_and_qc(file_path):
    """
    Parse file and apply quality control.
    Returns clean and rejected DataFrames.
    """

    df = parse_file(file_path)
    return run_quality_control(df)


   

    