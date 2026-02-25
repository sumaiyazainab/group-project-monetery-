```python
# =========================
# IMPORTS
# =========================

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
variant_analysis("DE_BSU_Pol_Batch_1.tsv", best_orf)
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
      <td>[R94G, L207L, E857G]</td>
      <td>[94, 207, 857]</td>
      <td>3</td>
      <td>1</td>
      <td>2</td>
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
      <td>[P83P, G260G, S340N, K464Q]</td>
      <td>[83, 260, 340, 464]</td>
      <td>4</td>
      <td>2</td>
      <td>2</td>
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
      <td>False</td>
    </tr>
  </tbody>
</table>
<p>301 rows × 14 columns</p>
</div>




```python
testing123 = variant_analysis("DE_BSU_Pol_Batch_1.tsv", best_orf)
```


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
      <td>[R94G, L207L, E857G]</td>
      <td>[94, 207, 857]</td>
      <td>3</td>
      <td>1</td>
      <td>2</td>
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
      <td>[P83P, G260G, S340N, K464Q]</td>
      <td>[83, 260, 340, 464]</td>
      <td>4</td>
      <td>2</td>
      <td>2</td>
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
      <td>False</td>
      <td>9.907855</td>
      <td>0.964086</td>
    </tr>
  </tbody>
</table>
<p>301 rows × 16 columns</p>
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
      <td>[R94G, L207L, E857G]</td>
      <td>[94, 207, 857]</td>
      <td>3</td>
      <td>1</td>
      <td>2</td>
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
      <td>[P83P, G260G, S340N, K464Q]</td>
      <td>[83, 260, 340, 464]</td>
      <td>4</td>
      <td>2</td>
      <td>2</td>
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
      <td>False</td>
      <td>9.907855</td>
      <td>0.922936</td>
    </tr>
  </tbody>
</table>
<p>301 rows × 16 columns</p>
</div>


