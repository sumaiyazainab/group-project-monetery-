#!/usr/bin/env python
# coding: utf-8

# ## TEJ'S CODE

# In[3]:


#Data Parsing & QC

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


# ## AMEERAH'S CODE

# In[4]:


import pandas as pd
import logging
from Bio.Seq import Seq

# =========================
# LOGGING
# =========================

logging.basicConfig(
    filename="phase6.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

# =========================
# QUALITY CONTROL
# =========================

# Dictionary for compulsory columns in uploaded data
ESSENTIAL_FIELDS = [
    "Plasmid_Variant_Index",
    "Parent_Plasmid_Variant",
    "Directed_Evolution_Generation",
    "Assembled_DNA_Sequence",
    "DNA_Quantification_fg",
    "Protein_Quantification_pg",
    "Control"
]


def parse_file(file_path):
    """
    Load TSV or JSON file.

    Separates I/O from biological logic for:
    - reproducibility
    - extensibility
    - format-agnosticism

    """
    if file_path.endswith(".tsv"):
        return pd.read_csv(file_path, sep="\t")
    elif file_path.endswith(".json"):
        return pd.read_json(file_path)
    else:
        raise ValueError("Unsupported file format")


def run_quality_control(df):
    """
    Validate required columns and remove incomplete records.

    Checks whether all columns listed in `ESSENTIAL_FIELDS`
    are present in the input DataFrame. If any required columns are missing,
    ValueError is raised. Then rows with missing values in any
    required columns are removed.

    Parameters
    ----------
    Input DataFrame to validate and clean.

    Returns
    -------
    tuple
        clean_df : DataFrame containing rows with all required fields present.
        rejected_df : DataFrame containing rows dropped due to missing
        values in required fields.

    Raises
    ------
    ValueError
        If one or more required columns defined in `ESSENTIAL_FIELDS`
        are not present in the DataFrame.

    """
    missing = [c for c in ESSENTIAL_FIELDS if c not in df.columns]
    if missing:
        raise ValueError(f"Missing essential columns: {missing}")

    clean_df = df.dropna(subset=ESSENTIAL_FIELDS)
    rejected_df = df.loc[~df.index.isin(clean_df.index)]
    return clean_df, rejected_df


def parse_and_qc(file_path):
    """
    Parse file and apply quality control checks.

    Parameters
    ----------
    file_path : str
        Path to the input file.

    Returns
    -------
    DataFrame
        Parsed and quality-controlled.

    """
    df = parse_file(file_path)
    return run_quality_control(df)


# =========================
# ORF FINDING BY SEED
# =========================

def best_approx_match_pos(seq, seed, max_mismatches=20):
    """
    Find best approximate match position of a seed sequence within a larger sequence.

    Slides a window across `seq` and computes the number of mismatches
    between each window and `seed`. Returns the position with the fewest
    mismatches, provided it does not exceed `max_mismatches`.

    Parameters
    ----------
    seq : str
        Full sequence to search within.
    seed : str
        Sub-sequence to match against.
    max_mismatches : int, optional
        Maximum mismatches allowed for a valid match (default is 20).

    Returns
    -------
    tuple
        (best_position, mismatch_count) if valid match is found;
        otherwise (None, mismatch_count).

    """
    L = len(seed)
    best_pos = None
    best_mm = L + 1

    for i in range(0, len(seq) - L + 1):
        window = seq[i:i+L]
        mm = sum(a != b for a, b in zip(window, seed))
        if mm < best_mm:
            best_mm = mm
            best_pos = i
            if mm == 0:
                return best_pos, 0

    if best_pos is not None and best_mm <= max_mismatches:
        return best_pos, best_mm
    return None, best_mm


def find_orf_by_seed(variant_seq, best_orf, seed_len=50, max_mismatches=20):
    """
    Identify open reading frame (ORF) in a circular plasmid sequence
    using a seed region from a reference (wild-type) ORF.

    Function searches both forward strand and reverse complement
    of the circularised variant sequence for best approximate match
    to the first `seed_len` nucleotides of the reference ORF. The match
    with the fewest mismatches (up to `max_mismatches`) is selected, and
    the full-length ORF sequence is extracted from that position.

    Parameters
    ----------
    variant_seq : str
        Nucleotide sequence of variant plasmid (assumed circular).
    best_orf : dict
        Dictionary containing reference ORF sequence under the key
        "orf_nt".
    seed_len : int, optional
        Length of 5' seed region from the reference ORF used for
        approximate matching (default is 50).
    max_mismatches : int, optional
        Maximum allowed mismatches when identifying the seed match
        (default is 20).

    Returns
    -------
    str
        Extracted nucleotide sequence of the identified ORF.

    Raises
    ------
    ValueError
        If the seed is not found on either strand, if invalid nucleotide
        characters are encountered, or if the extracted ORF length is
        not divisible by three.    

    """
    seq = variant_seq.upper()
    wt_seq = best_orf["orf_nt"].upper()
    seed = wt_seq[:seed_len]

    circular = seq + seq

    rev = circular[::-1]
    comp = {"A":"T","T":"A","G":"C","C":"G","N":"N"}
    try:
        circular_rev_comp = "".join(comp[b] for b in rev)
    except KeyError as e:
        raise ValueError(f"Invalid base in sequence: {e}")

    pos_f, mm_f = best_approx_match_pos(circular, seed, max_mismatches)
    pos_r, mm_r = best_approx_match_pos(circular_rev_comp, seed, max_mismatches)

    if pos_f is None and pos_r is None:
        raise ValueError("Seed not found on either strand")

    if pos_r is not None and (pos_f is None or mm_r <= mm_f):
        nt_seq = circular_rev_comp[pos_r:pos_r+len(wt_seq)]
    else:
        nt_seq = circular[pos_f:pos_f+len(wt_seq)]

    if len(nt_seq) % 3 != 0:
        raise ValueError("Extracted CDS not divisible by 3")

    return nt_seq


# =========================
# MUTATION CLASSIFICATION
# =========================

def classify_mutations(wt_cds, var_cds):
    """
    Compare two coding sequences (CDS) at codon resolution and classify mutations.

    Function translates both wild-type and variant CDS sequences and
    iterates codon-by-codon to determine whether each nucleotide change
    results in a synonymous (silent) or nonsynonymous (amino acid–altering)
    mutation. It also checks for premature stop codons in the variant
    protein sequence.

    Parameters
    ----------
    wt_cds : str
        Wild-type coding DNA sequence. Length must be divisible by 3.
    var_cds : str
        Variant coding DNA sequence aligned to the wild type and of equal length.

    Returns
    -------
    tuple
        total_mutations : int
            Total number of codon changes.
        synonymous_count : int
            Number of synonymous (silent) mutations.
        nonsynonymous_count : int
            Number of nonsynonymous (amino acid–changing) mutations.
        truncating : bool
            True if a premature stop codon is detected in the variant
            (excluding a terminal stop), otherwise False.

    Notes
    -----
    Function assumes both sequences are in-frame and properly aligned.

    """
    syn = 0
    nonsyn = 0
    trunc = False

    orf_aa = Seq(wt_cds).translate()
    var_prot = Seq(var_cds).translate()

    if "*" in var_prot[:-1]:
        trunc = True

    for i in range(0, len(wt_cds), 3):
        wt_codon = wt_cds[i:i+3]
        var_codon = var_cds[i:i+3]

        if wt_codon == var_codon:
            continue

        wt_aa = Seq(wt_codon).translate()
        var_aa = Seq(var_codon).translate()

        if wt_aa == var_aa:
            syn += 1
        else:
            nonsyn += 1

    return syn+nonsyn, syn, nonsyn, trunc



# =========================
# MUTATION LIST 
# =========================

def list_mutations(wt_cds, var_cds):
    """
    Generate amino acid–level mutation annotations between two CDS sequences.

    Function compares wild-type and variant sequences codon-by-codon,
    translates each codon, and records amino acid substitutions using
    standard notation (e.g., "W1L" for Trp→Leu at position 1).

    Parameters
    ----------
    wt_cds : str
        Wild-type coding DNA sequence. Length must be divisible by 3.
    var_cds : str
        Variant coding DNA sequence aligned to wild type and of equal length.

    Returns
    -------
    mutations : list of str
        List of amino acid substitution annotations in the format:
        "<WT_AA><position><VAR_AA>" (1-based indexing).
    aa_positions : list of int
        List of amino acid positions where mutations occur.

    Notes
    -----
    Only codons that differ at the nucleotide level are evaluated.
    Positions are reported using 1-based amino acid indexing.

    """
    muts = []
    aa_positions = []
    pos = 1

    for i in range(0, len(wt_cds), 3):
        wt_codon = wt_cds[i:i+3]
        var_codon = var_cds[i:i+3]

        if wt_codon != var_codon:
            wt_aa = str(Seq(wt_codon).translate())
            var_aa = str(Seq(var_codon).translate())
            muts.append(f"{wt_aa}{pos}{var_aa}")
            aa_positions.append(pos)

        pos += 1

    return muts, aa_positions


# =========================
# PHASE 6 PIPELINE
# =========================

def variant_analysis(data_path, best_orf):
    """
    Execute the Phase 6 variant analysis pipeline.

    Function processes plasmid sequence data by:
        1. Parsing and performing quality control on the input dataset.
        2. Identifying the wild-type ORF from the control plasmid.
        3. Extracting ORFs from each variant sequence using a seed-based search.
        4. Translating coding sequences to protein sequences.
        5. Classifying mutations relative to the wild type.
        6. Generating mutation annotations and summary statistics.
        7. Exporting results to CSV and JSON files.

    Parameters
    ----------
    data_path : str
        Path to the input dataset containing plasmid variant information.
    best_orf : dict
        Dictionary containing reference ORF information, including:
            "orf_nt" : wild-type nucleotide sequence
            "orf_aa" : wild-type amino acid sequence

    Returns
    -------
    DataFrame
            Containing the original input data along with:
            protein_seq : translated variant protein sequence
            mutations : list of amino acid substitutions
            aa_positions : mutated amino acid positions
            mutation_count : total number of mutations
            synonymous : count of synonymous mutations
            nonsynonymous : count of nonsynonymous mutations
            truncating : whether a premature stop codon is present

    Raises
    ------
    ValueError
        If wild-type ORF cannot be identified.

    Notes
    -----
    Function assumes plasmid sequences are circular and that
    all CDS sequences are in-frame and comparable to the wild type.

    Output files:
        - variant_analysis_output.csv
        - variant_analysis_output.json


    """
    df, rejected = parse_and_qc(data_path)

    wt_row = df[df["Control"] == True].iloc[0]

    wt_plasmid = wt_row["Assembled_DNA_Sequence"]

    wt_protein = best_orf["orf_aa"]

    wt_nt = find_orf_by_seed(wt_plasmid, best_orf)

    if isinstance(wt_nt, dict):
        raise ValueError(f"WT ORF could not be found: {wt_nt}")

    results = []

    for _, row in df.iterrows():

        logging.info(f"Processing variant {row['Plasmid_Variant_Index']}")

        var_nt = find_orf_by_seed(row["Assembled_DNA_Sequence"], best_orf)

        if isinstance(var_nt, dict):
            logging.warning(
                f"ORF not found for variant {row['Plasmid_Variant_Index']}: {var_nt}"
            )
            continue

        var_prot = str(Seq(var_nt).translate(to_stop=False))

        total, syn, nonsyn, trunc = classify_mutations(wt_nt, var_nt)

        mut_list, aa_positions = list_mutations(wt_nt, var_nt)

        results.append({
            "protein_seq": var_prot,
            "mutations": mut_list,
            "aa_positions": aa_positions,
            "mutation_count": total,
            "synonymous": syn,
            "nonsynonymous": nonsyn,
            "truncating": trunc
        })

    result_df = pd.concat([df.reset_index(drop=True),
                            pd.DataFrame(results)], axis=1)

    result_df.to_csv("variant_analysis_output.csv", index=False)

    result_df.to_json("variant_analysis_output.json", orient="records", indent=2)

    return result_df
# ===========================
# PHASE 7 — ACTIVITY SCORING
# ===========================

def identify_wt_rows(df):
    """
    Identify WT (baseline) control measurements using the Control column.

    Biological definition:
    ----------------------
    In prototype dataset, WT is not defined solely as generation 0.
    Instead, each generation has its own parental (baseline) construct.

    These baseline measurements are explicitly labelled:
        Control == True

    Therefore:
        WT rows = all rows where Control == True

    Returns
    -------
    DataFrame
        Subset of df containing only WT control rows.

    Raises
    ------
    ValueError
        If no control rows are found.

    """
    wt = df[df["Control"] == True]

    if len(wt) == 0:
        raise ValueError("No WT control rows detected (Control == True).")

    return wt


def compute_raw_activity(df):
    """
    Compute raw biochemical activity.

    Definition:
    -----------
        raw_activity = DNA_Quantification_fg / Protein_Quantification_pg

    Biological interpretation:
    --------------------------
    Estimates functional protein output per unit DNA template,
    correcting for variation in plasmid abundance or loading.

    This matches the experimental intent:
        activity ∝ protein produced per DNA input

    Returns
    -------
    DataFrame
        Copy of df with added column 'raw_activity'.

    """
    df = df.copy()

    df["raw_activity"] = (
        df["DNA_Quantification_fg"]/
        df["Protein_Quantification_pg"]
    )

    return df


def compute_global_wt_baseline(df):
    """
    Compute single global WT baseline activity.

    Definition:
    -----------
        global_baseline = mean(raw_activity) across ALL Control == True rows

    Biological meaning:
    -------------------
    This represents the average activity of all parental constructs
    measured across the experiment.

    This provides:
        • stable ancestral reference
        • common scale for cross-generation comparison

    Returns
    -------
    float
        Global baseline activity value.

    Raises
    ------
    ValueError
        If baseline is non-positive or undefined.

    """
    wt = identify_wt_rows(df)

    baseline = wt["raw_activity"].mean()

    if baseline <= 0:
        raise ValueError("Invalid global WT baseline (non-positive).")

    return baseline


def compute_generation_wt_baselines(df):
    """
    Compute WT baseline activity separately for each generation.

    Definition:
    -----------
        baseline[g] = mean(raw_activity) of rows where:
            Directed_Evolution_Generation == g
            AND Control == True

    Biological meaning:
    -------------------
    Each generation has its own parental construct.
    Mutants should be evaluated relative to their immediate parent.

    This enables:
        • correction for batch effects
        • within-generation fitness comparisons
        • local evolutionary step analysis

    Returns
    -------
    dictionary
        Mapping: generation -> baseline activity

    Raises
    ------
    ValueError
        If no generation baselines are computed.

    """
    wt = identify_wt_rows(df)

    baselines = (
        wt
        .groupby("Directed_Evolution_Generation")["raw_activity"]
        .mean()
        .to_dict()
    )

    if len(baselines) == 0:
        raise ValueError("No generation-specific WT baselines computed.")

    return baselines


def compute_activity_score(df, baseline_mode="global"):
    """
    Compute Activity Score using either global or generation-specific baseline.

    Parameters
    ----------
    DataFrame
        Must contain 'raw_activity' and generation metadata.

    baseline_mode : str
        Either:
            "global"     -> compare to global WT baseline
            "generation" -> compare to generation-specific WT baseline

    Formula:
    --------
        Activity_Score = raw_activity / baseline

    Where:
        baseline =
            global WT mean (if baseline_mode == "global")
            generation WT mean (if baseline_mode == "generation")

    Returns
    -------
    DataFrame
        Copy of df with added column 'Activity_Score'.

    Raises
    ------
    ValueError
        If baseline_mode is invalid.

    """
    df = df.copy()

    if baseline_mode == "global":
        baseline = compute_global_wt_baseline(df)
        df["Activity_Score"] = df["raw_activity"] / baseline

    elif baseline_mode == "generation":
        baselines = compute_generation_wt_baselines(df)

        df["Activity_Score"] = df.apply(
            lambda r: r["raw_activity"] /
            baselines[r["Directed_Evolution_Generation"]],
            axis=1
        )

    else:
        raise ValueError("baseline_mode must be 'global' or 'generation'.")

    return df


def activity_scoring(df, baseline_mode="global"):
    """
    Phase 7 — Activity Score Calculation

    Pipeline:
    ---------
    1. Compute raw biochemical activity
    2. Normalise by WT baseline
       (global or generation-specific)

    Inputs
    ------
    DataFrame
        Output of Phase 6 (mutation-annotated variants)

    baseline_mode : str
        "global" or "generation"

    Outputs
    -------
    DataFrame
        df with added columns:
            • raw_activity
            • Activity_Score

    Phase 7:
    -------------------------------
      + Identify WT control measurements
      + Compute baseline DNA & protein values
      + Define Activity Score formula
      + Normalise DNA yield by protein yield
      + Correct for WT baseline
      + Compute Activity Score for each variant
      + Store Activity Score in database

    """
    df = compute_raw_activity(df)

    df = compute_activity_score(df, baseline_mode)

    df.to_csv("activity_scores.csv", index=False)

    df.to_json("activity_scores.json", orient="records", indent=2)

    return df


# ## Ameerahs testing 

# In[17]:


best_orf = {'strand': '+',
'frame': 'f1',
'start_nt': 5070,
'end_nt': 7713,
'wraps_origin': False,
'length_nt': 2643,
'orf_nt': 'ATGACGGAACGTAAAAAACTGGTTCTGGTAGACGGCAACTCATTAGCCTATCGCGCTTTTTTTGCGTTGCCCCTGCTGTCTAATGACAAAGGAGTGCATACCAATGCTGTCTACGGCTTTGCGATGATTTTGATGAAAATGCTGGAAGATGAAAAACCCACCCACATGCTGGTTGCATTTGATGCAGGTAAAACGACGTTTCGTCATGGAACGTTTAAAGAATACAAAGGCGGTCGTCAGAAAACGCCGCCCGAACTGAGCGAACAGATGCCGTTTATCCGTGAACTGTTAGATGCATATCAAATCTCGCGCTATGAACTTGAACAATATGAGGCCGACGATATTATCGGTACCTTGGCAAAATCAGCGGAAAAAGATGGTTTTGAGGTTAAAGTATTTAGTGGGGACAAAGATCTTACACAACTGGCGACCGATAAGACCACCGTAGCCATCACCCGCAAAGGTATTACCGATGTAGAGTTCTATACCCCGGAACACGTCAAAGAAAAATACGGTTTAACTCCGGAACAGATTATTGATATGAAAGGCCTGATGGGCGACTCTAGCGACAATATCCCGGGGGTGCCGGGTGTTGGTGAAAAAACTGCCATTAAGCTGTTAAAACAGTTTGATTCCGTAGAAAAACTGCTGGAATCAATTGATGAGGTGTCTGGAAAGAAATTGAAAGAGAAATTAGAAGAATTCAAGGACCAGGCCTTAATGTCTAAGGAACTTGCAACAATCATGACCGATGCGCCGATTGAAGTCAGTGTTAGCGGGCTGGAGTACCAGGGGTTTAATCGCGAACAAGTGATCGCAATTTTTAAAGATCTGGGCTTTAACACGCTGCTTGAACGCCTTGGGGAAGATAGCGCCGAAGCCGAACAGGATCAGTCTTTGGAGGACATTAACGTAAAGACCGTTACTGATGTTACTTCAGACATTTTAGTCTCCCCGTCTGCCTTCGTGGTAGAGCAGATTGGCGATAATTATCACGAAGAACCGATCTTGGGTTTCAGCATCGTGAACGAAACTGGCGCGTACTTTATTCCGAAAGATATTGCGGTGGAAAGCGAGGTTTTCAAAGAATGGGTAGAGAATGACGAGCAAAAGAAATGGGTGTTTGATTCAAAACGCGCGGTGGTGGCGCTGCGTTGGCAGGGCATTGAGCTGAAAGGGGCGGAATTCGATACCCTGCTTGCAGCGTACATTATTAACCCGGGGAACAGCTACGATGATGTGGCGTCGGTAGCCAAAGATTATGGCTTACACATTGTTAGCTCTGACGAATCCGTGTATGGCAAAGGCGCCAAGCGCGCGGTGCCTAGCGAAGATGTGTTATCAGAACACCTTGGTCGCAAAGCCCTGGCCATTCAGAGTCTGCGTGAAAAGCTGGTGCAAGAGCTTGAGAACAATGATCAGTTGGAACTGTTTGAAGAGTTGGAGATGCCGCTGGCGCTGATCCTGGGAGAAATGGAAAGCACGGGTGTGAAAGTGGATGTCGACCGCTTAAAACGCATGGGCGAAGAGTTAGGAGCCAAGCTGAAAGAGTATGAGGAGAAAATTCACGAAATTGCTGGCGAACCGTTTAACATCAACTCCCCTAAGCAGTTGGGCGTGATTCTGTTTGAAAAAATCGGTCTTCCGGTGGTGAAGAAAACTAAAACCGGGTATTCCACCAGTGCCGATGTGCTTGAGAAACTGGCGGATAAACATGATATCGTAGATTACATCTTGCAGTACCGCCAGATTGGAAAACTGCAGTCTACATATATTGAAGGCTTACTTAAAGTGACGCGCCCGGACAGTCACAAGGTCCACACCCGCTTTAATCAGGCACTGACCCAGACTGGCCGTCTTAGCAGCACCGATCCGAATCTGCAGAACATTCCGATTCGCTTGGAAGAGGGTCGCAAAATTCGTCAGGCCTTTGTGCCGAGCGAAAAGGACTGGTTAATCTTCGCGGCCGATTATAGCCAGATTGAACTGCGTGTACTGGCCCATATTAGTAAAGATGAGAATCTGATCGAGGCGTTTACCAACGATATGGACATTCATACTAAAACCGCTATGGATGTTTTTCATGTAGCCAAAGATGAAGTCACTTCTGCGATGCGCCGTCAGGCTAAAGCCGTTAACTTTGGCATTGTATACGGGATTAGCGACTACGGCCTTAGCCAGAACCTGGGCATTACGCGTAAAGAAGCGGGGGCGTTTATCGATCGTTACCTGGAATCGTTCCAGGGCGTGAAAGCATACATGGAGGACAGTGTGCAGGAAGCGAAACAAAAAGGCTACGTTACCACCTTGATGCACCGTCGTCGCTATATCCCAGAACTGACGAGTCGCAATTTTAATATTCGCAGCTTCGCCGAGCGCACGGCGATGAACACACCGATCCAGGGGTCTGCGGCAGATATTATTAAAAAAGCAATGATCGACATGGCCGCAAAGCTTAAAGAAAAGCAACTGAAAGCGCGCTTATTACTTCAGGTTCATGATGAGCTGATCTTTGAAGCGCCAAAAGAGGAAATTGAAATCCTGGAAAAACTGGTACCAGAAGTAATGGAACATGCGCTTGCGTTGGATGTGCCGTTAAAAGTTGACTTCGCCTCAGGACCGAGCTGGTATGATGCCAAATAA',
'length_aa': 880,
'orf_aa': 'MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMKMLEDEKPTHMLVAFDAGKTTFRHGTFKEYKGGRQKTPPELSEQMPFIRELLDAYQISRYELEQYEADDIIGTLAKSAEKDGFEVKVFSGDKDLTQLATDKTTVAITRKGITDVEFYTPEHVKEKYGLTPEQIIDMKGLMGDSSDNIPGVPGVGEKTAIKLLKQFDSVEKLLESIDEVSGKKLKEKLEEFKDQALMSKELATIMTDAPIEVSVSGLEYQGFNREQVIAIFKDLGFNTLLERLGEDSAEAEQDQSLEDINVKTVTDVTSDILVSPSAFVVEQIGDNYHEEPILGFSIVNETGAYFIPKDIAVESEVFKEWVENDEQKKWVFDSKRAVVALRWQGIELKGAEFDTLLAAYIINPGNSYDDVASVAKDYGLHIVSSDESVYGKGAKRAVPSEDVLSEHLGRKALAIQSLREKLVQELENNDQLELFEELEMPLALILGEMESTGVKVDVDRLKRMGEELGAKLKEYEEKIHEIAGEPFNINSPKQLGVILFEKIGLPVVKKTKTGYSTSADVLEKLADKHDIVDYILQYRQIGKLQSTYIEGLLKVTRPDSHKVHTRFNQALTQTGRLSSTDPNLQNIPIRLEEGRKIRQAFVPSEKDWLIFAADYSQIELRVLAHISKDENLIEAFTNDMDIHTKTAMDVFHVAKDEVTSAMRRQAKAVNFGIVYGISDYGLSQNLGITRKEAGAFIDRYLESFQGVKAYMEDSVQEAKQKGYVTTLMHRRRYIPELTSRNFNIRSFAERTAMNTPIQGSAADIIKKAMIDMAAKLKEKQLKARLLLQVHDELIFEAPKEEIEILEKLVPEVMEHALALDVPLKVDFASGPSWYDAK',
'alignment_score': 4444.0,
'score_fraction': 1.0,
'filtered_orf_count': 10,
'total_orf_count': 190,
'status': 'match_found'}


# In[4]:


testing123 = variant_analysis("DE_BSU_Pol_Batch_1.tsv", best_orf)
testing123


# In[5]:


#testing123 = variant_analysis("DE_BSU_Pol_Batch_1.tsv", best_orf)
#activity_scoring(testing123, baseline_mode = "global")


# In[6]:


mutation_count_df = activity_scoring(testing123, baseline_mode = "generation")
mutation_count_df


# ## VISUALISATION - SUMAIYAS CODE

# In[5]:


from __future__ import annotations

from typing import List, Dict, Tuple, Optional

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import plotly.graph_objects as go

from sklearn.decomposition import PCA
from sklearn.feature_extraction.text import CountVectorizer
from sklearn.neighbors import KNeighborsRegressor

import ast


# In[6]:


# 1) Helpers: required columns + validation


ESSENTIAL_FIELDS = [
    "Plasmid_Variant_Index",
    "Parent_Plasmid_Variant",
    "Directed_Evolution_Generation",
    "DNA_Quantification_fg",
    "Protein_Quantification_pg",
    "Control",
]

MUTATION_FIELDS_EXPECTED = [
    "mutation_count",   # produced by Phase 6 variant_analysis
    "synonymous",
    "nonsynonymous",
    "truncating",
]

ACTIVITY_FIELDS_EXPECTED = [
    "raw_activity",     # produced by Phase 7 activity_scoring
    "Activity_Score",
]


def _require_columns(df: pd.DataFrame, required: List[str], name: str = "DataFrame") -> None:
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"{name} is missing required columns: {missing}")


def ensure_mutation_count(df: pd.DataFrame) -> pd.DataFrame:
    """
    Ensure df contains a mutation_count column.

    Expected:
      - If Phase 6 ran, df already has 'mutation_count'
      - Otherwise we add mutation_count=NaN (so tables still render)

    Returns:
      Copy of df with mutation_count present.
    """
    out = df.copy()
    if "mutation_count" not in out.columns:
        out["mutation_count"] = np.nan
    return out


# In[7]:


def ensure_mutation_count(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()

    # Keep if already exists
    if "mutation_count" in out.columns:
        out["mutation_count"] = pd.to_numeric(out["mutation_count"], errors="coerce")
        return out

    # Try derive from mutations list
    if "mutations" in out.columns:
        def _to_list(x):
            if isinstance(x, list):
                return x
            if isinstance(x, str):
                try:
                    return ast.literal_eval(x)  # handles "['A12V']"
                except Exception:
                    return None
            return None

        muts = out["mutations"].apply(_to_list)
        out["mutation_count"] = muts.apply(lambda x: len(x) if isinstance(x, list) else np.nan)
        return out

    # Try derive from aa_positions list
    if "aa_positions" in out.columns:
        def _to_list(x):
            if isinstance(x, list):
                return x
            if isinstance(x, str):
                try:
                    return ast.literal_eval(x)
                except Exception:
                    return None
            return None

        pos = out["aa_positions"].apply(_to_list)
        out["mutation_count"] = pos.apply(lambda x: len(x) if isinstance(x, list) else np.nan)
        return out

    out["mutation_count"] = np.nan
    return out


# In[11]:


def top_performers_table(
    df: pd.DataFrame,
    n: int = 10,
    exclude_controls: bool = True,
    essential_fields: Optional[List[str]] = None,
    score_col: str = "Activity_Score",
) -> pd.DataFrame:
    """
    Return the top N performing variants by Activity Score.

    Output columns:
      essential fields + Activity_Score +
      mutations/aa_positions (if present) +
      mutation summary columns (mutation_count, synonymous, nonsynonymous, truncating if present)
    """
    if essential_fields is None:
        essential_fields = ESSENTIAL_FIELDS

    _require_columns(df, [score_col], name="top_performers_table input")
    _require_columns(
        df,
        [c for c in essential_fields if c != "Control"],
        name="top_performers_table input",
    )

    out = ensure_mutation_count(df)

    if exclude_controls and "Control" in out.columns:
        out = out[out["Control"] == False].copy()

    out = out.sort_values(score_col, ascending=False).copy()

    # Include mutation detail columns if available
    mutation_detail_cols = ["mutations", "aa_positions"]
    mutation_detail_cols = [c for c in mutation_detail_cols if c in out.columns]

    # Include mutation summary columns if available
    mutation_summary_cols = ["mutation_count", "synonymous", "nonsynonymous", "truncating"]
    mutation_summary_cols = [c for c in mutation_summary_cols if c in out.columns]

    cols = (
        [c for c in essential_fields if c in out.columns]
        + [score_col]
        + mutation_detail_cols
        + mutation_summary_cols
    )

    top = out[cols].head(n).reset_index(drop=True)

    # tidy formatting
    top[score_col] = (
        pd.to_numeric(top[score_col], errors="coerce")
        .astype(float)
        .round(4)
    )

    return top


# In[12]:


# 3) ACTIVITY DISTRIBUTION PLOT (by generation)

def plot_activity_distribution_by_generation(
    df: pd.DataFrame,
    score_col: str = "Activity_Score",
    exclude_controls: bool = True,
    kind: str = "violin",  # "violin" or "box"
    figsize: Tuple[int, int] = (10, 5),
    show_means: bool = True,
    show_medians: bool = True,
    save_path: Optional[str] = None,   
    dpi: int = 300    
) -> None:
    """
    Plot Activity Score distribution per generation.

    Matplotlib-only so it's easy to embed as PNG in Flask templates.
    """
    _require_columns(df, ["Directed_Evolution_Generation", score_col], name="plot_activity_distribution_by_generation input")

    plot_df = df.copy()
    if exclude_controls and "Control" in plot_df.columns:
        plot_df = plot_df[plot_df["Control"] == False].copy()

    gens = sorted(plot_df["Directed_Evolution_Generation"].dropna().unique().tolist())
    data = [
        plot_df.loc[plot_df["Directed_Evolution_Generation"] == g, score_col].dropna().values
        for g in gens
    ]

    plt.figure(figsize=figsize)

    if kind.lower() == "violin":
        plt.violinplot(data, showmeans=show_means, showmedians=show_medians)
    elif kind.lower() == "box":
        plt.boxplot(data, showmeans=show_means)
    else:
        raise ValueError("kind must be 'violin' or 'box'.")

    plt.xticks(np.arange(1, len(gens) + 1), [str(g) for g in gens])
    plt.xlabel("Generation")
    plt.ylabel("WT-normalised Activity Score")
    plt.title("Activity Score Distribution by Generation")
    plt.grid(True, alpha=0.3)
    plt.tight_layout()

    if save_path:
        if not save_path.endswith(".png"):
            save_path += ".png"
        plt.savefig(save_path, dpi=dpi, bbox_inches="tight")

    plt.show()


# In[15]:


# 4) 3D ACTIVITY LANDSCAPE (Plotly)

def build_2d_sequence_embedding(
    df: pd.DataFrame,
    seq_col: str = "protein_seq",
    method: str = "pca",
    kmer_k: int = 3,
    n_components: int = 2,
    random_state: int = 42,
) -> pd.DataFrame:
    """
    Convert sequences into a 2D coordinate system.

    Default aligns with your old approach:
      - k-mer count vectors (CountVectorizer char n-grams)
      - PCA -> Dim1, Dim2

    Returns:
      Copy of df with added columns: Dim1, Dim2
    """
    _require_columns(df, [seq_col], name="build_2d_sequence_embedding input")

    out = df.copy()
    out[seq_col] = out[seq_col].astype(str)

    vectorizer = CountVectorizer(analyzer="char", ngram_range=(kmer_k, kmer_k))
    X = vectorizer.fit_transform(out[seq_col])

    if method.lower() == "pca":
        reducer = PCA(n_components=n_components, random_state=random_state)
        coords = reducer.fit_transform(X.toarray())
    else:
        raise ValueError("Only method='pca' is implemented here (fast + deterministic).")

    out["Dim1"] = coords[:, 0]
    out["Dim2"] = coords[:, 1]
    return out


def make_3d_activity_landscape(
    df: pd.DataFrame,
    score_col: str = "Activity_Score",
    seq_col: str = "protein_seq",
    exclude_controls: bool = True,
    kmer_k: int = 3,
    grid_size: int = 70,
    knn_neighbors: int = 5,
    knn_weights: str = "distance",
    title: str = "3D Activity Landscape (Sequence Diversity → PCA → Activity Score)",   
    save_html: Optional[str] = None, 
) -> go.Figure:
    """
    Create a 3D topography-like surface where:
      - X, Y = 2D sequence embedding (PCA over k-mers)
      - Z = Activity Score
      - Surface = KNN-regressed activity over grid
      - Points = actual variants

    Returns:
      plotly.graph_objects.Figure (embed into Flask with fig.to_html())
    """
    _require_columns(df, [score_col, seq_col], name="make_3d_activity_landscape input")

    plot_df = df.copy()
    if exclude_controls and "Control" in plot_df.columns:
        plot_df = plot_df[plot_df["Control"] == False].copy()

    plot_df = plot_df.dropna(subset=[seq_col, score_col]).copy()

    # Build Dim1/Dim2 from sequences
    plot_df = build_2d_sequence_embedding(
        plot_df,
        seq_col=seq_col,
        method="pca",
        kmer_k=kmer_k,
        n_components=2,
        random_state=42,
    )

    # Training data
    xy = plot_df[["Dim1", "Dim2"]].values
    z = plot_df[score_col].astype(float).values

    # Grid
    x_lin = np.linspace(plot_df["Dim1"].min(), plot_df["Dim1"].max(), grid_size)
    y_lin = np.linspace(plot_df["Dim2"].min(), plot_df["Dim2"].max(), grid_size)
    xx, yy = np.meshgrid(x_lin, y_lin)
    grid_xy = np.c_[xx.ravel(), yy.ravel()]

    # KNN regression surface
    knn = KNeighborsRegressor(n_neighbors=knn_neighbors, weights=knn_weights)
    knn.fit(xy, z)
    zz = knn.predict(grid_xy).reshape(grid_size, grid_size)

    fig = go.Figure()

    # Surface
    fig.add_trace(go.Surface(
        x=xx, y=yy, z=zz,
        surfacecolor=zz,
        opacity=0.85,
        colorbar=dict(title="Activity Score")
    ))

    # Points
    hover_text = []
    for _, r in plot_df.iterrows():
        vid = r.get("Plasmid_Variant_Index", "")
        gen = r.get("Directed_Evolution_Generation", "")
        score = r.get(score_col, np.nan)
        mut = r.get("mutation_count", np.nan)
        hover_text.append(
            f"Variant {vid}<br>Gen {gen}<br>Score {float(score):.3f}<br>Mutations {mut}"
        )

    fig.add_trace(go.Scatter3d(
        x=plot_df["Dim1"], y=plot_df["Dim2"], z=plot_df[score_col],
        mode="markers",
        marker=dict(size=3, opacity=0.7),
        text=hover_text,
        hoverinfo="text",
        name="Variants"
    ))

    fig.update_layout(
        title=title,
        scene=dict(
            xaxis_title="Sequence diversity (PCA Dim 1)",
            yaxis_title="Sequence diversity (PCA Dim 2)",
            zaxis_title="Activity Score"
        ),
        width=900,
        height=650,
        margin=dict(l=0, r=0, b=0, t=50)
    )

    if save_html:
        if not save_html.endswith(".html"):
            save_html += ".html"

        fig.write_html(
            save_html,
            include_plotlyjs="cdn",   
            full_html=True
        )


    return fig


# ## sumaiyas testing

# In[18]:


testing123 = variant_analysis("DE_BSU_Pol_Batch_1.tsv", best_orf)
mutation_count_df = activity_scoring(testing123, baseline_mode = "generation")
mutation_count_df


# In[19]:


qc_df = ensure_mutation_count(mutation_count_df)


# In[20]:


top10 = top_performers_table(qc_df, n=10)
top10


# In[21]:


plot_activity_distribution_by_generation(
    mutation_count_df,
    score_col="Activity_Score",
    exclude_controls=True,
    kind="violin",      # or "box"
    figsize=(10, 5)
)


# In[22]:


fig = make_3d_activity_landscape(
    mutation_count_df,
    score_col="Activity_Score",
    seq_col="protein_seq",
    exclude_controls=True,
    kmer_k=3
)
fig.show()


# In[ ]:




