import pandas as pd

# Essential fields
ESSENTIAL_FIELDS = [
    "Plasmid_Variant_Index",
    "Parent_Plasmid_Variant",
    "Directed_Evolution_Generation",
    "Assembled_DNA_Sequence",
    "DNA_Quantification_fg",
    "Protein_Quantification_pg"
]


def parse_file(file_path):
    """Load TSV or JSON file"""
    if file_path.endswith(".tsv"):
        df = pd.read_csv(file_path, sep="\t")
    elif file_path.endswith(".json"):
        df = pd.read_json(file_path)
    else:
        raise ValueError("Unsupported file format")

    return df


def run_quality_control(df):
    """Validate schema + remove incomplete rows"""

    missing_columns = [col for col in ESSENTIAL_FIELDS if col not in df.columns]
    if missing_columns:
        raise ValueError(f"Missing essential columns: {missing_columns}")

    clean_df = df.dropna(subset=ESSENTIAL_FIELDS)
    rejected_df = df.loc[~df.index.isin(clean_df.index)]

    return clean_df, rejected_df


def parse_and_qc(file_path):
    df = parse_file(file_path)
    return run_quality_control(df)
