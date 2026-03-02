import pandas as pd
import logging
from Bio.Seq import Seq



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

def variant_analysis(df, best_orf):
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


    return df
