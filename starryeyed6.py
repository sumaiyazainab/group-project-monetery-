#!/usr/bin/env python
# coding: utf-8

# In[3]:


# =========================
# IMPORTS
# =========================

import pandas as pd
import logging
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner

aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 2
aligner.mismatch_score = -1
aligner.open_gap_score = -5
aligner.extend_gap_score = -1


# In[2]:


# =========================
# LOGGING
# =========================

logging.basicConfig(
    filename="variantanalysis.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)


# In[3]:


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


# In[4]:


# ===========================
# PHASE 6 - VARIANT ANALYSIS
# ===========================


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


def extract_cds_to_first_stop(circular_seq: str, start_pos: int, max_cds_nt: int = 20000) -> str:
    """
    Extract coding sequence (CDS) from circular DNA sequence beginning
    at a specified start position and ending at the first in-frame stop codon.

    Stop codon is included in the returned CDS.

    Parameters
    ----------
    circular_seq : str
        DNA sequence (must already be circularised, e.g., seq + seq).
    start_pos : int
        Index of first nucleotide of the start codon.
    max_cds_nt : int, optional
        Maximum number of nucleotides to scan (default 20000).

    Returns
    -------
    str
        Nucleotide sequence from start codon to first in-frame stop (inclusive).

    Raises
    ------
    ValueError
        If no in-frame stop codon is found within max_cds_nt.

    """
    STOP_CODONS = {"TAA", "TAG", "TGA"}

    # Ensure we have enough sequence to walk; circular_seq should already be doubled.
    seq = circular_seq[start_pos:start_pos + max_cds_nt]

    # Walk codons; include the stop codon in the returned CDS
    for i in range(0, len(seq) - 2, 3):
        codon = seq[i:i+3]
        if codon in STOP_CODONS:
            cds = seq[:i+3]
            return cds

    raise ValueError("No in-frame stop codon found within max_cds_nt")



def find_orf_by_seed_and_stop(variant_seq, best_orf, seed_len=50, max_mismatches=20):
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
        start_pos = pos_r
        use_seq = circular_rev_comp
    else:
        start_pos = pos_f
        use_seq = circular

    cds = extract_cds_to_first_stop(use_seq, start_pos)

    if len(cds) % 3 != 0:
        raise ValueError("Extracted CDS not divisible by 3")

    return cds


def translate_with_alt_start_with_stop(sequence: str, start_codons: set[str] | None = None) -> str:
    """
    Translate a nucleotide coding sequence into a protein sequence.

    The translation:
    - Includes stop codons ('*') in output.
    - Forces first amino acid to Methionine ('M') if the first codon
      is in the allowed start_codons set (e.g., ATG, GTG, TTG).

    Parameters
    ----------
    sequence : str
        DNA coding sequence (must be in-frame).
    start_codons : set[str] or None, optional
        Codons treated as valid alternative start codons.
        Default: {"ATG", "GTG", "TTG"}.

    Returns
    -------
    str
        Translated amino acid sequence including stop codons.

    """

    seq = sequence.upper()

    # if user did not provide start_codons selection, use default start_codons
    if start_codons is None:
        start_codons = {"ATG", "GTG", "TTG"}

    # using biopython translate module, translate sequence up until the stop codon - do not include *
    aa_seq = str(Seq(seq).translate(to_stop=False))
    first_codon = seq[:3]

    # if the first_codon in aa_seq is in the selected start_codons, force the first codon into M for UniProt reference sequence alignment
    if first_codon in start_codons and aa_seq:
        aa_seq = "M" + aa_seq[1:]

    return aa_seq




def global_align_proteins(wt_prot: str, var_prot: str):
    """
    Perform global protein alignment between wild-type and variant sequences
    using Biopython's PairwiseAligner.

    Scoring scheme:
        match = +2
        mismatch = -1
        gap open = -5
        gap extend = -1

    Parameters
    ----------
    wt_prot : str
        Wild-type amino acid sequence.
    var_prot : str
        Variant amino acid sequence.

    Returns
    -------
    tuple[str, str]
        (aligned_wt_sequence, aligned_variant_sequence)

    """

    alignments = aligner.align(wt_prot, var_prot)
    best_alignment = alignments[0]

    # Build gapped sequences manually
    wt_aln = []
    var_aln = []

    wt_index = 0
    var_index = 0

    for (wt_start, wt_end), (var_start, var_end) in zip(
        best_alignment.aligned[0],
        best_alignment.aligned[1],
    ):
        # Add gaps before this aligned block
        while wt_index < wt_start:
            wt_aln.append(wt_prot[wt_index])
            var_aln.append("-")
            wt_index += 1

        while var_index < var_start:
            wt_aln.append("-")
            var_aln.append(var_prot[var_index])
            var_index += 1

        # Add aligned block
        while wt_index < wt_end and var_index < var_end:
            wt_aln.append(wt_prot[wt_index])
            var_aln.append(var_prot[var_index])
            wt_index += 1
            var_index += 1

    # Add trailing gaps
    while wt_index < len(wt_prot):
        wt_aln.append(wt_prot[wt_index])
        var_aln.append("-")
        wt_index += 1

    while var_index < len(var_prot):
        wt_aln.append("-")
        var_aln.append(var_prot[var_index])
        var_index += 1

    return "".join(wt_aln), "".join(var_aln)


def build_codon_aligned_cds(wt_aln: str, var_aln: str, wt_nt_orf: str, var_nt_orf: str):
    """
    Convert gapped protein alignment (wt_aln/var_aln) into codon-aligned CDS strings.
    Protein gap '-' becomes codon gap '---'.

    """
    wt_nt_orf = wt_nt_orf.upper()
    var_nt_orf = var_nt_orf.upper()

    if len(wt_aln) != len(var_aln):
        raise ValueError("Protein alignments must have equal length")

    if len(wt_nt_orf) % 3 != 0 or len(var_nt_orf) % 3 != 0:
        raise ValueError("ORF nucleotide sequences must be divisible by 3")

    wt_i = 0  # codon index in WT ORF
    var_i = 0 # codon index in VAR ORF

    wt_out = []
    var_out = []

    for a, b in zip(wt_aln, var_aln):
        # WT side
        if a == "-":
            wt_out.append("---")
        else:
            wt_out.append(wt_nt_orf[wt_i*3: wt_i*3+3])
            wt_i += 1

        # VAR side
        if b == "-":
            var_out.append("---")
        else:
            var_out.append(var_nt_orf[var_i*3: var_i*3+3])
            var_i += 1

    wt_cds_gap = "".join(wt_out)
    var_cds_gap = "".join(var_out)

    if len(wt_cds_gap) != len(var_cds_gap):
        raise ValueError("Internal error: gapped CDS lengths not equal (should never happen)")

    return wt_cds_gap, var_cds_gap


def classify_mutations_codon_aligned(wt_cds: str, var_cds: str):
    """
    Classify mutations between codon-aligned wild-type and variant CDS sequences.

    Determines:
        - synonymous substitutions
        - nonsynonymous substitutions
        - amino acid insertions
        - amino acid deletions
        - internal stop codons (truncation)

    Parameters
    ----------
    wt_cds : str
        Codon-aligned wild-type nucleotide sequence (may contain '---').
    var_cds : str
        Codon-aligned variant nucleotide sequence (may contain '---').

    Returns
    -------
    tuple
        (
            mutation_count,
            synonymous_count,
            nonsynonymous_count,
            truncation_flag,
            insertions_count,
            deletions_count
        )

    """

    wt_cds = wt_cds.upper()
    var_cds = var_cds.upper()

    if len(wt_cds) != len(var_cds):
        raise ValueError("Codon-aligned CDS must have equal length")
    if len(wt_cds) % 3 != 0:
        raise ValueError("Aligned CDS length must be divisible by 3")

    syn = 0
    nonsyn = 0
    insertions_aa = 0
    deletions_aa = 0
    trunc = False

    # truncation
    # var_aln comes from previous step before alignment
    var_nt_nogap = var_cds.replace("-", "")
    if len(var_nt_nogap) % 3 != 0:
        raise ValueError("Ungapped variant CDS not divisible by 3 (frameshift or bad gapping)")

    var_prot = str(Seq(var_nt_nogap).translate(to_stop=False))


    # Ignore terminal stop(s), flag only internal stops
    var_prot_no_terminal_stop = var_prot.rstrip("*")
    trunc = ("*" in var_prot_no_terminal_stop)

    for i in range(0, len(wt_cds), 3):
        wt_codon = wt_cds[i:i+3]
        var_codon = var_cds[i:i+3]

        if wt_codon == var_codon:
            continue
        if wt_codon == "---" and var_codon != "---":
            insertions_aa += 1
            continue
        if wt_codon != "---" and var_codon == "---":
            deletions_aa += 1
            continue
        if wt_codon != var_codon:
            wt_aa = str(Seq(wt_codon).translate())
            var_aa = str(Seq(var_codon).translate())
            if wt_aa == var_aa:
                syn += 1
            else:
                nonsyn += 1

    mutation_count = syn + nonsyn + insertions_aa + deletions_aa
    return mutation_count, syn, nonsyn, trunc, insertions_aa, deletions_aa


def extract_mutation_annotations(wt_aln: str, var_aln: str):
    """
    Generate amino acid mutation annotations from aligned protein sequences.

    Returns:
        mutations: list of mutation strings (e.g. ["T599S"])
        positions: list of 1-based amino acid positions

    """
    mutations = []
    positions = []

    wt_pos = 0  # counts WT residues (ignore gaps)

    for a, b in zip(wt_aln, var_aln):
        if a != "-":
            wt_pos += 1

        # substitution
        if a != "-" and b != "-" and a != b:
            mutations.append(f"{a}{wt_pos}{b}")
            positions.append(wt_pos)

        # deletion
        elif a != "-" and b == "-":
            mutations.append(f"{a}{wt_pos}del")
            positions.append(wt_pos)

        # insertion
        elif a == "-" and b != "-":
            mutations.append(f"ins{wt_pos}{b}")
            positions.append(wt_pos)

    return mutations, positions


# =========================
# PHASE 6 PIPELINE - updated
# =========================

def variant_analysis(data_path, best_orf, seed_len=80, max_mismatches=8):
    """
    Execute Phase 6 variant analysis pipeline.

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
            Contains original input data along with:
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

    # identity wt control is in row 1 - control parent plasmid
    wt_row = df[df["Control"] == True].iloc[0]
    wt_plasmid = wt_row["Assembled_DNA_Sequence"]

    # extract the wt orf
    wt_nt = find_orf_by_seed_and_stop(wt_plasmid, best_orf, seed_len=seed_len, max_mismatches=max_mismatches)
    if isinstance(wt_nt, dict):
        raise ValueError(f"WT ORF could not be found: {wt_nt}")

    # translate wt_nt from extracted orf
    wt_aa_seq = translate_with_alt_start_with_stop(wt_nt)

    # find the variant sequence orf
    results = []
    for _, row in df.iterrows():
        logging.info(f"Processing variant {row['Plasmid_Variant_Index']}")
        try:
            var_plasmid = row["Assembled_DNA_Sequence"]
            var_nt_orf = find_orf_by_seed_and_stop(var_plasmid, best_orf, seed_len=seed_len, max_mismatches=max_mismatches)

            # translate
            var_aa_seq = translate_with_alt_start_with_stop(var_nt_orf)

            # protein alignment (THIS WAS MISSING)
            wt_aln, var_aln = global_align_proteins(wt_aa_seq, var_aa_seq)

            # extract mutation annotations
            mutations, aa_positions = extract_mutation_annotations(wt_aln, var_aln)

            # build codon-aligned CDS
            wt_cds_gap, var_cds_gap = build_codon_aligned_cds(
                wt_aln,
                var_aln,
                wt_nt,
                var_nt_orf
            )

            # classify mutations (including possible indels)
            mutation_count, syn, nonsyn, trunc, insertions_aa, deletions_aa = classify_mutations_codon_aligned(wt_cds_gap, var_cds_gap)

            results.append({
                "protein_seq": var_aa_seq,
                "mutations": mutations,
                "aa_positions": aa_positions,
                "mutation_count": mutation_count,
                "synonymous": syn,
                "nonsynonymous": nonsyn,
                "insertions": insertions_aa,
                "deletions": deletions_aa,
                "truncating": trunc
            })

        except Exception as e:
            logging.warning(
                f"Variant {row['Plasmid_Variant_Index']} failed: {e}"
            )
            results.append({
                "protein_seq": None,
                "mutations": None,
                "aa_positions": None,
                "mutation_count": None,
                "synonymous": None,
                "nonsynonymous": None,
                "insertions": None,
                "deletions": None,
                "truncating": None,
                "error": str(e)
            })

    # merge results
    result_df = pd.concat([df.reset_index(drop=True),
                            pd.DataFrame(results).reset_index(drop=True)], axis=1)

    result_df.to_csv("variant_analysis_output.csv", index=False)
    result_df.to_json("variant_analysis_output.json", orient="records", indent=2)

    return result_df


# In[5]:


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

    Output files:
        - activity_scores.csv
        - activity_scores.json

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


# In[6]:


# =========================
# TESTING
# =========================


# In[7]:


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


# In[8]:


variant_analysis("DE_BSU_Pol_Batch_1.tsv", best_orf, seed_len=80, max_mismatches=8)


# In[9]:


testing123 = variant_analysis("DE_BSU_Pol_Batch_1.tsv", best_orf, seed_len=80, max_mismatches=8)


# In[12]:


testing123.tail(30)


# In[10]:


activity_scoring(testing123, baseline_mode = "global")


# In[11]:


activity_scoring(testing123, baseline_mode = "generation")

