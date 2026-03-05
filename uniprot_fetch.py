import pandas as pd
import requests

def fetch_uniprot_information(accession_id: str) -> dict:
    
    """
    Fetch protein information from the UniProt database.
    This function takes an UniProt Accession ID input and fetches the associated protein information from the UniProt database as a JSON report for downstream data extraction and analysis.
    :param accession_id: UniProt accession ID (e.g. P05067)
    :type accession_id: str
    :return: Parsed UniProt protein information as a JSON dictionary
    :rtype: dict
    """

    # get JSON for the accession id UniProt page with a timeout enforced where code fails fast if UniProt does not respond
    try:
        response = requests.get(f"https://rest.uniprot.org/uniprotkb/{accession_id}.json", timeout = 10)
        if response.status_code == 200:
            # successful request: return parsed JSON for downstream analysis
            return response.json()
        # accession ID is valid but not found in UniProt
        elif response.status_code == 404:
            raise ValueError(f"The UniProt accession '{accession_id}' does not exist in the UniProt database. Please check your input and try again.")
        # malformed accession ID (e.g., invalid format)
        elif response.status_code == 400:
            raise ValueError("The UniProt ID format you entered is invalid. Please check your input and try again.")
        # UniProt server-side failure
        elif 500 <= response.status_code < 600:
            raise RuntimeError("We could not retrieve the protein due to a UniProt server error. Please try again later.")
        else:
            # unexpected HTTP responses
            raise RuntimeError("An unexpected error has occured while retrieving the protein.")
    # handles UniProt request timeouts to avoid blocking downstream pipelines
    except requests.exceptions.Timeout:
        raise RuntimeError("We could not connect to UniProt in time. Please try again.")
    # handles all other network-related issues (e.g. connection error)
    except requests.exceptions.RequestException:
        raise RuntimeError("An unexpected signal/network error occurred. Please try again later.")
    
    

def extract_uniprot_aa_sequence(uniprot_record: dict) -> str:

    """
    Extract amino acid sequence from a UniProt protein information record.
    This function retrieves the protein amino acid sequence from a UniProt JSON record from the UniProt REST API. A RuntimeError is raised if the sequence is missing.
    :param uniprot_record: Parsed UniProt protein information as a JSON dictionary
    :type uniprot_record: dict
    :return: Amino acid sequence of the protein of interest
    :rtype: str
    """

    # extract the amino acid sequence from JSON under sequence, value
    try:
        uniprot_protein_seq = uniprot_record["sequence"]["value"]
        return uniprot_protein_seq
    # handles malformed or incomplete UniProt records where the protein sequence field is missing or has an invalid structure
    except (KeyError, TypeError):
        raise RuntimeError("The UniProt record is missing a protein sequence.")
    

def extract_uniprot_features(uniprot_record: dict) -> pd.DataFrame:

    """
    Extract location and function of all features from a UniProt protein information record.
    This function extracts UniProt feature annotations for a target protein and flattens into a DataFrame for downstream analysis and SQL storage.
    :param uniprot_record: Parsed UniProt protein information as a JSON dictionary
    :type uniprot_record: dict
    :return: Dataframe of protein features (including feature type, description, start_pos and end_pos)
    :rtype: pd.DataFrame
    """

    features = []
    # for each feature, if no feature return empty list, if there is an annotated feature, get the start and end positions from the location block and append the details into features list
    for x in uniprot_record.get("features", []):
        location = x.get("location", {})
        start = location.get("start",{}).get("value")
        end = location.get("end",{}).get("value")
        features.append({
            "type":  x.get("type"),
            "description": x.get("description"),
            "start_pos": start,
            "end_pos": end
            })
    features_df = pd.DataFrame(features, columns=["type", "description", "start_pos", "end_pos"])
    return features_df


def extract_key_domains(features_df: pd.DataFrame) -> pd.DataFrame:

    """
    Filters the location and function of key functional domains (Domain features only) of a protein from a DataFrame of UniProt annotated features.
    This function filters UniProt key domain annotations from a base DataFrame with all annotated features for downstream analysis.
    :param features_df: Dataframe of protein features (including feature type, description, start_pos and end_pos)
    :type features_df: pd.DataFrame
    :return: Dataframe of key functional protein domains (including description, start_pos and end_pos)
    :rtype: pd.DataFrame
    """

    # filter domain-only features into a new key domains dataframe - case-insensitive matching is applied to ensure robust feature selection
    # prevent errors if the "type" column contains null values
    key_domains_df = features_df[features_df["type"].fillna("").str.upper() == "DOMAIN"]
    return key_domains_df
