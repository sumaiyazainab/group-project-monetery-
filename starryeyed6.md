```python
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
```


```python
# =========================
# LOGGING
# =========================

logging.basicConfig(
    filename="variantanalysis.log",
    level=logging.INFO,
    format="%(asctime)s - %(levelname)s - %(message)s"
)

```


```python
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
```


```python
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
```


```python
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
```


```python
# =========================
# TESTING
# =========================
```


```python
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
```


```python
variant_analysis("DE_BSU_Pol_Batch_1.tsv", best_orf, seed_len=80, max_mismatches=8)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Plasmid_Variant_Index</th>
      <th>Parent_Plasmid_Variant</th>
      <th>Directed_Evolution_Generation</th>
      <th>Assembled_DNA_Sequence</th>
      <th>DNA_Quantification_fg</th>
      <th>Protein_Quantification_pg</th>
      <th>Control</th>
      <th>protein_seq</th>
      <th>mutations</th>
      <th>aa_positions</th>
      <th>mutation_count</th>
      <th>synonymous</th>
      <th>nonsynonymous</th>
      <th>insertions</th>
      <th>deletions</th>
      <th>truncating</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>-1</td>
      <td>0</td>
      <td>CTGGTTGCATTTGATGCAGGTAAAACGACGTTTCGTCATGGAACGT...</td>
      <td>521.01</td>
      <td>52.98</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>CTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGC...</td>
      <td>744.46</td>
      <td>61.06</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[T599S]</td>
      <td>[599]</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>0</td>
      <td>1</td>
      <td>GTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCTTTA...</td>
      <td>735.16</td>
      <td>56.29</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[T769P]</td>
      <td>[769]</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3</td>
      <td>0</td>
      <td>1</td>
      <td>GACAAGCTGTGACCGTCTCCGGGAGCTGCATGTGTCAGAGGTTTTC...</td>
      <td>729.84</td>
      <td>57.86</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[R94G, E857G]</td>
      <td>[94, 857]</td>
      <td>3</td>
      <td>1</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4</td>
      <td>0</td>
      <td>1</td>
      <td>CATCGGTCGAGATCCCGGTGCCTAATGAGTGAGCTAACTTACATTA...</td>
      <td>725.82</td>
      <td>55.71</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[S340N, K464Q]</td>
      <td>[340, 464]</td>
      <td>4</td>
      <td>2</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>296</th>
      <td>296</td>
      <td>0</td>
      <td>10</td>
      <td>ATCCGTGTATGGCAAAGGCGCCAAGCGCGCGGTGCCTAGCGAAGAT...</td>
      <td>426.45</td>
      <td>42.03</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>297</th>
      <td>297</td>
      <td>0</td>
      <td>10</td>
      <td>ATATAAATCAGCATCCATGTTGGAATTTAATCGCGGCCTAGAGCAA...</td>
      <td>782.60</td>
      <td>68.86</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>298</th>
      <td>298</td>
      <td>0</td>
      <td>10</td>
      <td>CTGTTTGAAAAAATCGGTCTTCCGGTGGTGAAGAAAACTAAAACCG...</td>
      <td>733.14</td>
      <td>60.81</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>299</th>
      <td>299</td>
      <td>0</td>
      <td>10</td>
      <td>CTGATAAAGCGGGCCATGTTAAGGGCGGTTTTTTCCTGTTTGGTCA...</td>
      <td>551.63</td>
      <td>54.08</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>300</th>
      <td>300</td>
      <td>0</td>
      <td>10</td>
      <td>GTCCAGCTCGTTGAGTTTCTCCAGAAGCGTTAATGTCTGGCTTCTG...</td>
      <td>524.72</td>
      <td>52.96</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
<p>301 rows × 16 columns</p>
</div>




```python
testing123 = variant_analysis("DE_BSU_Pol_Batch_1.tsv", best_orf, seed_len=80, max_mismatches=8)
```


```python
testing123.tail(30)
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Plasmid_Variant_Index</th>
      <th>Parent_Plasmid_Variant</th>
      <th>Directed_Evolution_Generation</th>
      <th>Assembled_DNA_Sequence</th>
      <th>DNA_Quantification_fg</th>
      <th>Protein_Quantification_pg</th>
      <th>Control</th>
      <th>protein_seq</th>
      <th>mutations</th>
      <th>aa_positions</th>
      <th>mutation_count</th>
      <th>synonymous</th>
      <th>nonsynonymous</th>
      <th>insertions</th>
      <th>deletions</th>
      <th>truncating</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>271</th>
      <td>271</td>
      <td>247</td>
      <td>10</td>
      <td>CCGCCAACACCCGCTGACGCGCCCTGACGGGCTTGTCTGCTCCCGG...</td>
      <td>2000.00</td>
      <td>52.58</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, E51A, E106D, I115V, V256del, S291G, Q32...</td>
      <td>[31, 51, 106, 115, 256, 291, 326, 397, 405, 44...</td>
      <td>21</td>
      <td>2</td>
      <td>15</td>
      <td>0</td>
      <td>4</td>
      <td>False</td>
    </tr>
    <tr>
      <th>272</th>
      <td>272</td>
      <td>244</td>
      <td>10</td>
      <td>GTTACTGATGTTACTTCAGACATTTTAGTCTCCCCGTCTGCCTTCG...</td>
      <td>2000.00</td>
      <td>51.51</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, E51A, T71S, I115V, V256del, S291G, Q326...</td>
      <td>[31, 51, 71, 115, 256, 291, 326, 356, 397, 427...</td>
      <td>22</td>
      <td>4</td>
      <td>15</td>
      <td>0</td>
      <td>3</td>
      <td>False</td>
    </tr>
    <tr>
      <th>273</th>
      <td>273</td>
      <td>242</td>
      <td>10</td>
      <td>ATCCTGCGATGCAGATCCGGAACATAATGGTGCAGGGCGCTGACTT...</td>
      <td>2000.00</td>
      <td>44.59</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, E51A, T71S, I115V, T144A, V256del, S291...</td>
      <td>[31, 51, 71, 115, 144, 256, 291, 326, 397, 422...</td>
      <td>24</td>
      <td>5</td>
      <td>16</td>
      <td>0</td>
      <td>3</td>
      <td>False</td>
    </tr>
    <tr>
      <th>274</th>
      <td>274</td>
      <td>257</td>
      <td>10</td>
      <td>TTATCGGTACCTTGGCAACATCAGCGGAAAAAGATGGTTTTGAGGT...</td>
      <td>2000.00</td>
      <td>57.50</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, E51A, I115V, K121T, V256del, S291G, Q32...</td>
      <td>[31, 51, 115, 121, 256, 291, 326, 397, 443, 45...</td>
      <td>30</td>
      <td>3</td>
      <td>17</td>
      <td>7</td>
      <td>3</td>
      <td>False</td>
    </tr>
    <tr>
      <th>275</th>
      <td>275</td>
      <td>257</td>
      <td>10</td>
      <td>CATGTTAAGGGCGGTTTTTTCCTGTTTGGTCACTGATGCCTCCGTG...</td>
      <td>2000.00</td>
      <td>55.80</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, E51A, I115V, K121T, V256del, S291G, Q32...</td>
      <td>[31, 51, 115, 121, 256, 291, 326, 397, 443, 45...</td>
      <td>29</td>
      <td>3</td>
      <td>16</td>
      <td>7</td>
      <td>3</td>
      <td>False</td>
    </tr>
    <tr>
      <th>276</th>
      <td>276</td>
      <td>244</td>
      <td>10</td>
      <td>CGGTAAGCGGCAGGGTCGGAACAGGAGAGCGCACGAGGGAGCTTCC...</td>
      <td>2000.00</td>
      <td>57.97</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, E51A, T71S, I115V, V256del, S291G, Q326...</td>
      <td>[31, 51, 71, 115, 256, 291, 326, 397, 427, 443...</td>
      <td>22</td>
      <td>4</td>
      <td>15</td>
      <td>0</td>
      <td>3</td>
      <td>False</td>
    </tr>
    <tr>
      <th>277</th>
      <td>277</td>
      <td>244</td>
      <td>10</td>
      <td>TAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTCCTTCTAGT...</td>
      <td>2000.00</td>
      <td>53.71</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, E51A, T71S, I115V, K121N, V256del, S291...</td>
      <td>[31, 51, 71, 115, 121, 256, 291, 326, 397, 427...</td>
      <td>22</td>
      <td>4</td>
      <td>15</td>
      <td>0</td>
      <td>3</td>
      <td>False</td>
    </tr>
    <tr>
      <th>278</th>
      <td>278</td>
      <td>254</td>
      <td>10</td>
      <td>ATTCACCACCCTGAATTGACTCTCTTCCGGGCGCTATCATGCCATA...</td>
      <td>2000.00</td>
      <td>59.36</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLASRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[Y17S, G31V, I115V, V256del, S291G, Q326L, D39...</td>
      <td>[17, 31, 115, 256, 291, 326, 397, 443, 497, 56...</td>
      <td>19</td>
      <td>1</td>
      <td>16</td>
      <td>0</td>
      <td>2</td>
      <td>False</td>
    </tr>
    <tr>
      <th>279</th>
      <td>279</td>
      <td>245</td>
      <td>10</td>
      <td>AGCATCGCAGTGGGAACGATGCCCTCATTCAGCATTTGCATGGTTT...</td>
      <td>2000.00</td>
      <td>61.98</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, E51A, T71S, I115V, V256del, S291G, Q326...</td>
      <td>[31, 51, 71, 115, 256, 291, 326, 397, 427, 443...</td>
      <td>22</td>
      <td>4</td>
      <td>15</td>
      <td>0</td>
      <td>3</td>
      <td>False</td>
    </tr>
    <tr>
      <th>280</th>
      <td>280</td>
      <td>247</td>
      <td>10</td>
      <td>GACGCGCCCTGTAGCGGCGCATTAAGCGCGGCGGGTGTGGTGGTTA...</td>
      <td>2000.00</td>
      <td>58.74</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, E51A, E106D, I115V, V256del, S291G, Q32...</td>
      <td>[31, 51, 106, 115, 256, 291, 326, 397, 443, 49...</td>
      <td>21</td>
      <td>3</td>
      <td>14</td>
      <td>0</td>
      <td>4</td>
      <td>False</td>
    </tr>
    <tr>
      <th>281</th>
      <td>281</td>
      <td>244</td>
      <td>10</td>
      <td>GCTTCCCATACAATCGATAGATTGTCGCACCTGATTGCCCGACATT...</td>
      <td>2000.00</td>
      <td>55.20</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, E51A, T71S, I115V, V256del, S291G, Q326...</td>
      <td>[31, 51, 71, 115, 256, 291, 326, 397, 427, 443...</td>
      <td>23</td>
      <td>4</td>
      <td>16</td>
      <td>0</td>
      <td>3</td>
      <td>False</td>
    </tr>
    <tr>
      <th>282</th>
      <td>282</td>
      <td>257</td>
      <td>10</td>
      <td>CCCTTTAGGGTTCCGATTTAGTGCTTTACGGCACCTCGACCCCAAA...</td>
      <td>2000.00</td>
      <td>53.11</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, E51A, I115V, K121T, V256del, S291G, Q32...</td>
      <td>[31, 51, 115, 121, 256, 291, 326, 397, 443, 45...</td>
      <td>28</td>
      <td>3</td>
      <td>15</td>
      <td>7</td>
      <td>3</td>
      <td>False</td>
    </tr>
    <tr>
      <th>283</th>
      <td>283</td>
      <td>249</td>
      <td>10</td>
      <td>ATTTTCACCTGAATCAGGATATTCTTCTAATACCTGGAATGCTGTT...</td>
      <td>1992.55</td>
      <td>64.49</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLASRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[Y17S, G31V, I115V, V256del, S291G, Q326L, F37...</td>
      <td>[17, 31, 115, 256, 291, 326, 375, 397, 443, 49...</td>
      <td>19</td>
      <td>2</td>
      <td>15</td>
      <td>0</td>
      <td>2</td>
      <td>False</td>
    </tr>
    <tr>
      <th>284</th>
      <td>284</td>
      <td>260</td>
      <td>10</td>
      <td>CGGTAATGGCGCGCATTGCGCCCAGCGCCATCTGATCGTTGGCAAC...</td>
      <td>1980.07</td>
      <td>51.02</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, E51A, T71S, I115V, K226Q, V256del, S291...</td>
      <td>[31, 51, 71, 115, 226, 256, 291, 326, 397, 427...</td>
      <td>22</td>
      <td>3</td>
      <td>16</td>
      <td>0</td>
      <td>3</td>
      <td>False</td>
    </tr>
    <tr>
      <th>285</th>
      <td>285</td>
      <td>251</td>
      <td>10</td>
      <td>ATGGCGCGCATTGCGCCCAGCGCCATCTGATCGTTGGCAACCAGCA...</td>
      <td>1979.58</td>
      <td>63.16</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLASRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[Y17S, G31V, I115V, V256del, S291G, Q326L, D39...</td>
      <td>[17, 31, 115, 256, 291, 326, 397, 443, 497, 52...</td>
      <td>18</td>
      <td>2</td>
      <td>14</td>
      <td>0</td>
      <td>2</td>
      <td>False</td>
    </tr>
    <tr>
      <th>286</th>
      <td>286</td>
      <td>242</td>
      <td>10</td>
      <td>ATGACGAGCAAAAGAAATGGGTGTTTGATTCAAAACGCGCGGTGGT...</td>
      <td>1978.82</td>
      <td>59.78</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, E51A, T71S, I115V, T144A, V256del, S291...</td>
      <td>[31, 51, 71, 115, 144, 256, 291, 326, 397, 427...</td>
      <td>22</td>
      <td>4</td>
      <td>15</td>
      <td>0</td>
      <td>3</td>
      <td>False</td>
    </tr>
    <tr>
      <th>287</th>
      <td>287</td>
      <td>259</td>
      <td>10</td>
      <td>GCGTACATTATTAACCCGGGGAACAGCTACGATGATGTGGCGTCGG...</td>
      <td>1974.48</td>
      <td>59.89</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, I115V, D126Y, L216P, V256del, S291G, D3...</td>
      <td>[31, 115, 126, 216, 256, 291, 302, 326, 397, 4...</td>
      <td>18</td>
      <td>1</td>
      <td>15</td>
      <td>0</td>
      <td>2</td>
      <td>False</td>
    </tr>
    <tr>
      <th>288</th>
      <td>288</td>
      <td>250</td>
      <td>10</td>
      <td>CCCTGATAGACGGTTTTTCGCCCTTTGACGTTGGAGTCCACGTTCT...</td>
      <td>1968.83</td>
      <td>48.11</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRPFFALPLLSNDKVVLTNAVYGFAMILMK...</td>
      <td>[A19P, G31V, H33L, E51A, I115V, V256del, S291G...</td>
      <td>[19, 31, 33, 51, 115, 256, 291, 326, 397, 443,...</td>
      <td>35</td>
      <td>4</td>
      <td>20</td>
      <td>7</td>
      <td>4</td>
      <td>False</td>
    </tr>
    <tr>
      <th>289</th>
      <td>289</td>
      <td>260</td>
      <td>10</td>
      <td>GAATTGTGAGCGGATAACAATTCCCCTCTAGAAATAATTTTGTTTA...</td>
      <td>1962.97</td>
      <td>58.67</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, E51A, T71S, I115V, V256del, S291G, K306...</td>
      <td>[31, 51, 71, 115, 256, 291, 306, 326, 397, 427...</td>
      <td>22</td>
      <td>3</td>
      <td>16</td>
      <td>0</td>
      <td>3</td>
      <td>False</td>
    </tr>
    <tr>
      <th>290</th>
      <td>290</td>
      <td>259</td>
      <td>10</td>
      <td>CATGACGGAACGTAAAAAACTGGTTCTGGTAGACGGCAACTCATTA...</td>
      <td>1961.51</td>
      <td>56.20</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, I115V, L216P, V256del, S291G, D302Y, Q3...</td>
      <td>[31, 115, 216, 256, 291, 302, 326, 397, 443, 4...</td>
      <td>18</td>
      <td>2</td>
      <td>14</td>
      <td>0</td>
      <td>2</td>
      <td>False</td>
    </tr>
    <tr>
      <th>291</th>
      <td>291</td>
      <td>254</td>
      <td>10</td>
      <td>CCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTGCTGAAAGG...</td>
      <td>1955.88</td>
      <td>56.61</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLASRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[Y17S, G31V, I115V, A143T, V256del, S291G, Q32...</td>
      <td>[17, 31, 115, 143, 256, 291, 326, 397, 443, 49...</td>
      <td>19</td>
      <td>1</td>
      <td>16</td>
      <td>0</td>
      <td>2</td>
      <td>False</td>
    </tr>
    <tr>
      <th>292</th>
      <td>292</td>
      <td>254</td>
      <td>10</td>
      <td>TCGATAGATTGTCGCACCTGATTGCCCGACATTATCGCGAGCCCAT...</td>
      <td>1955.16</td>
      <td>65.62</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLASRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[Y17S, G31V, I115V, V256del, S291G, Q326L, D39...</td>
      <td>[17, 31, 115, 256, 291, 326, 397, 425, 443, 49...</td>
      <td>19</td>
      <td>1</td>
      <td>16</td>
      <td>0</td>
      <td>2</td>
      <td>False</td>
    </tr>
    <tr>
      <th>293</th>
      <td>293</td>
      <td>253</td>
      <td>10</td>
      <td>AGAAAAATACGGTTTAACTCCGGAACAGATTATTGATATGAAAGGC...</td>
      <td>1946.45</td>
      <td>49.69</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLASRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[Y17S, G31V, I115V, V256del, V258A, S291G, Q32...</td>
      <td>[17, 31, 115, 256, 258, 291, 326, 397, 443, 49...</td>
      <td>20</td>
      <td>1</td>
      <td>17</td>
      <td>0</td>
      <td>2</td>
      <td>False</td>
    </tr>
    <tr>
      <th>294</th>
      <td>294</td>
      <td>243</td>
      <td>10</td>
      <td>GTAACTGGCTTCAGCAGAGCGCAGATACCAAATACTGTCCTTCTAG...</td>
      <td>1940.33</td>
      <td>52.82</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[G31V, E51A, T71S, I115V, V256del, S291G, Q326...</td>
      <td>[31, 51, 71, 115, 256, 291, 326, 367, 397, 427...</td>
      <td>22</td>
      <td>3</td>
      <td>16</td>
      <td>0</td>
      <td>3</td>
      <td>False</td>
    </tr>
    <tr>
      <th>295</th>
      <td>295</td>
      <td>241</td>
      <td>10</td>
      <td>GCTTGCAGCGTACATTATTAACCCGGGGAACAGCTACGATGATGTG...</td>
      <td>1934.29</td>
      <td>48.66</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLASRAFFALPLLSNDKVVHTNAVYGFAMILMK...</td>
      <td>[Y17S, G31V, I115V, V256del, I272L, D290A, S29...</td>
      <td>[17, 31, 115, 256, 272, 290, 291, 326, 360, 39...</td>
      <td>22</td>
      <td>2</td>
      <td>18</td>
      <td>0</td>
      <td>2</td>
      <td>False</td>
    </tr>
    <tr>
      <th>296</th>
      <td>296</td>
      <td>0</td>
      <td>10</td>
      <td>ATCCGTGTATGGCAAAGGCGCCAAGCGCGCGGTGCCTAGCGAAGAT...</td>
      <td>426.45</td>
      <td>42.03</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>297</th>
      <td>297</td>
      <td>0</td>
      <td>10</td>
      <td>ATATAAATCAGCATCCATGTTGGAATTTAATCGCGGCCTAGAGCAA...</td>
      <td>782.60</td>
      <td>68.86</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>298</th>
      <td>298</td>
      <td>0</td>
      <td>10</td>
      <td>CTGTTTGAAAAAATCGGTCTTCCGGTGGTGAAGAAAACTAAAACCG...</td>
      <td>733.14</td>
      <td>60.81</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>299</th>
      <td>299</td>
      <td>0</td>
      <td>10</td>
      <td>CTGATAAAGCGGGCCATGTTAAGGGCGGTTTTTTCCTGTTTGGTCA...</td>
      <td>551.63</td>
      <td>54.08</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
    <tr>
      <th>300</th>
      <td>300</td>
      <td>0</td>
      <td>10</td>
      <td>GTCCAGCTCGTTGAGTTTCTCCAGAAGCGTTAATGTCTGGCTTCTG...</td>
      <td>524.72</td>
      <td>52.96</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
    </tr>
  </tbody>
</table>
</div>




```python
activity_scoring(testing123, baseline_mode = "global")
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Plasmid_Variant_Index</th>
      <th>Parent_Plasmid_Variant</th>
      <th>Directed_Evolution_Generation</th>
      <th>Assembled_DNA_Sequence</th>
      <th>DNA_Quantification_fg</th>
      <th>Protein_Quantification_pg</th>
      <th>Control</th>
      <th>protein_seq</th>
      <th>mutations</th>
      <th>aa_positions</th>
      <th>mutation_count</th>
      <th>synonymous</th>
      <th>nonsynonymous</th>
      <th>insertions</th>
      <th>deletions</th>
      <th>truncating</th>
      <th>raw_activity</th>
      <th>Activity_Score</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>-1</td>
      <td>0</td>
      <td>CTGGTTGCATTTGATGCAGGTAAAACGACGTTTCGTCATGGAACGT...</td>
      <td>521.01</td>
      <td>52.98</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>9.834088</td>
      <td>0.956908</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>CTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGC...</td>
      <td>744.46</td>
      <td>61.06</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[T599S]</td>
      <td>[599]</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>12.192270</td>
      <td>1.186371</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>0</td>
      <td>1</td>
      <td>GTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCTTTA...</td>
      <td>735.16</td>
      <td>56.29</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[T769P]</td>
      <td>[769]</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>13.060224</td>
      <td>1.270827</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3</td>
      <td>0</td>
      <td>1</td>
      <td>GACAAGCTGTGACCGTCTCCGGGAGCTGCATGTGTCAGAGGTTTTC...</td>
      <td>729.84</td>
      <td>57.86</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[R94G, E857G]</td>
      <td>[94, 857]</td>
      <td>3</td>
      <td>1</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>12.613896</td>
      <td>1.227397</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4</td>
      <td>0</td>
      <td>1</td>
      <td>CATCGGTCGAGATCCCGGTGCCTAATGAGTGAGCTAACTTACATTA...</td>
      <td>725.82</td>
      <td>55.71</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[S340N, K464Q]</td>
      <td>[340, 464]</td>
      <td>4</td>
      <td>2</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>13.028541</td>
      <td>1.267744</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>296</th>
      <td>296</td>
      <td>0</td>
      <td>10</td>
      <td>ATCCGTGTATGGCAAAGGCGCCAAGCGCGCGGTGCCTAGCGAAGAT...</td>
      <td>426.45</td>
      <td>42.03</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>10.146324</td>
      <td>0.987290</td>
    </tr>
    <tr>
      <th>297</th>
      <td>297</td>
      <td>0</td>
      <td>10</td>
      <td>ATATAAATCAGCATCCATGTTGGAATTTAATCGCGGCCTAGAGCAA...</td>
      <td>782.60</td>
      <td>68.86</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>11.365089</td>
      <td>1.105882</td>
    </tr>
    <tr>
      <th>298</th>
      <td>298</td>
      <td>0</td>
      <td>10</td>
      <td>CTGTTTGAAAAAATCGGTCTTCCGGTGGTGAAGAAAACTAAAACCG...</td>
      <td>733.14</td>
      <td>60.81</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>12.056241</td>
      <td>1.173135</td>
    </tr>
    <tr>
      <th>299</th>
      <td>299</td>
      <td>0</td>
      <td>10</td>
      <td>CTGATAAAGCGGGCCATGTTAAGGGCGGTTTTTTCCTGTTTGGTCA...</td>
      <td>551.63</td>
      <td>54.08</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>10.200259</td>
      <td>0.992538</td>
    </tr>
    <tr>
      <th>300</th>
      <td>300</td>
      <td>0</td>
      <td>10</td>
      <td>GTCCAGCTCGTTGAGTTTCTCCAGAAGCGTTAATGTCTGGCTTCTG...</td>
      <td>524.72</td>
      <td>52.96</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>9.907855</td>
      <td>0.964086</td>
    </tr>
  </tbody>
</table>
<p>301 rows × 18 columns</p>
</div>




```python
activity_scoring(testing123, baseline_mode = "generation")
```




<div>
<style scoped>
    .dataframe tbody tr th:only-of-type {
        vertical-align: middle;
    }

    .dataframe tbody tr th {
        vertical-align: top;
    }

    .dataframe thead th {
        text-align: right;
    }
</style>
<table border="1" class="dataframe">
  <thead>
    <tr style="text-align: right;">
      <th></th>
      <th>Plasmid_Variant_Index</th>
      <th>Parent_Plasmid_Variant</th>
      <th>Directed_Evolution_Generation</th>
      <th>Assembled_DNA_Sequence</th>
      <th>DNA_Quantification_fg</th>
      <th>Protein_Quantification_pg</th>
      <th>Control</th>
      <th>protein_seq</th>
      <th>mutations</th>
      <th>aa_positions</th>
      <th>mutation_count</th>
      <th>synonymous</th>
      <th>nonsynonymous</th>
      <th>insertions</th>
      <th>deletions</th>
      <th>truncating</th>
      <th>raw_activity</th>
      <th>Activity_Score</th>
    </tr>
  </thead>
  <tbody>
    <tr>
      <th>0</th>
      <td>0</td>
      <td>-1</td>
      <td>0</td>
      <td>CTGGTTGCATTTGATGCAGGTAAAACGACGTTTCGTCATGGAACGT...</td>
      <td>521.01</td>
      <td>52.98</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>9.834088</td>
      <td>1.000000</td>
    </tr>
    <tr>
      <th>1</th>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>CTGACTTGAGCGTCGATTTTTGTGATGCTCGTCAGGGGGGCGGAGC...</td>
      <td>744.46</td>
      <td>61.06</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[T599S]</td>
      <td>[599]</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>12.192270</td>
      <td>1.261237</td>
    </tr>
    <tr>
      <th>2</th>
      <td>2</td>
      <td>0</td>
      <td>1</td>
      <td>GTTCGCCGGCTTTCCCCGTCAAGCTCTAAATCGGGGGCTCCCTTTA...</td>
      <td>735.16</td>
      <td>56.29</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[T769P]</td>
      <td>[769]</td>
      <td>1</td>
      <td>0</td>
      <td>1</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>13.060224</td>
      <td>1.351023</td>
    </tr>
    <tr>
      <th>3</th>
      <td>3</td>
      <td>0</td>
      <td>1</td>
      <td>GACAAGCTGTGACCGTCTCCGGGAGCTGCATGTGTCAGAGGTTTTC...</td>
      <td>729.84</td>
      <td>57.86</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[R94G, E857G]</td>
      <td>[94, 857]</td>
      <td>3</td>
      <td>1</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>12.613896</td>
      <td>1.304853</td>
    </tr>
    <tr>
      <th>4</th>
      <td>4</td>
      <td>0</td>
      <td>1</td>
      <td>CATCGGTCGAGATCCCGGTGCCTAATGAGTGAGCTAACTTACATTA...</td>
      <td>725.82</td>
      <td>55.71</td>
      <td>False</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[S340N, K464Q]</td>
      <td>[340, 464]</td>
      <td>4</td>
      <td>2</td>
      <td>2</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>13.028541</td>
      <td>1.347746</td>
    </tr>
    <tr>
      <th>...</th>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
      <td>...</td>
    </tr>
    <tr>
      <th>296</th>
      <td>296</td>
      <td>0</td>
      <td>10</td>
      <td>ATCCGTGTATGGCAAAGGCGCCAAGCGCGCGGTGCCTAGCGAAGAT...</td>
      <td>426.45</td>
      <td>42.03</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>10.146324</td>
      <td>0.945149</td>
    </tr>
    <tr>
      <th>297</th>
      <td>297</td>
      <td>0</td>
      <td>10</td>
      <td>ATATAAATCAGCATCCATGTTGGAATTTAATCGCGGCCTAGAGCAA...</td>
      <td>782.60</td>
      <td>68.86</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>11.365089</td>
      <td>1.058680</td>
    </tr>
    <tr>
      <th>298</th>
      <td>298</td>
      <td>0</td>
      <td>10</td>
      <td>CTGTTTGAAAAAATCGGTCTTCCGGTGGTGAAGAAAACTAAAACCG...</td>
      <td>733.14</td>
      <td>60.81</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>12.056241</td>
      <td>1.123062</td>
    </tr>
    <tr>
      <th>299</th>
      <td>299</td>
      <td>0</td>
      <td>10</td>
      <td>CTGATAAAGCGGGCCATGTTAAGGGCGGTTTTTTCCTGTTTGGTCA...</td>
      <td>551.63</td>
      <td>54.08</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>10.200259</td>
      <td>0.950174</td>
    </tr>
    <tr>
      <th>300</th>
      <td>300</td>
      <td>0</td>
      <td>10</td>
      <td>GTCCAGCTCGTTGAGTTTCTCCAGAAGCGTTAATGTCTGGCTTCTG...</td>
      <td>524.72</td>
      <td>52.96</td>
      <td>True</td>
      <td>MTERKKLVLVDGNSLAYRAFFALPLLSNDKGVHTNAVYGFAMILMK...</td>
      <td>[]</td>
      <td>[]</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>0</td>
      <td>False</td>
      <td>9.907855</td>
      <td>0.922936</td>
    </tr>
  </tbody>
</table>
<p>301 rows × 18 columns</p>
</div>


