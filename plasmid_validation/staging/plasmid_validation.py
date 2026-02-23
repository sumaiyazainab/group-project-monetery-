
from Bio import Align
from Bio.Align import substitution_matrices

def length_compatibility(orf_aa_len: int, uniprot_len: int) -> bool:

    """
    Check amino acid length compatibility between candidate ORFs and UniProt sequence.
    This function evaluates whether a candidate orf sequence length is within ±20% of the uniprot sequence length, and is used as
    a fast pre-filtering step to remove too short or long ORFs before computationally expensive alignment.
    :param orf_aa_len: Amino acid length of the candidate ORF
    :type orf_aa_len: int
    :param uniprot_len: Amino acid length of the UniProt reference protein
    :type uniprot_len: int
    :return: True if ORF length is within ±20% of the UniProt length, otherwise False
    :rtype: bool
    """
    
    # prevents comparison if UniProt length is invalid
    if uniprot_len <= 0:
        return False
    
    # calculate the corresponding min and max length threshold for the uniprot length
    min_length = int(0.8 * uniprot_len)
    max_length = int(1.2 * uniprot_len)

    # return if ORF length is within ±20% of the UniProt length
    return min_length <= orf_aa_len <= max_length

def filter_orf_length(orfs: list[dict], uniprot_protein_seq: str) -> list[dict]:
    
    """
    Filter list of candidates depending on amino acid length compatibility with UniProt sequence.
    This function filters for candidate ORF sequences by length within ±20% of the reference uniprot. ORFs with missing 
    amino acid length information are skipped.
    :param orfs: List of dictionaries of candidate ORFs produced by candidate ORF detection
    :type orfs: list[dict]
    :param uniprot_protein_seq: UniProt reference protein amino acid sequence
    :type uniprot_protein_seq: str
    :return: Filtered list of ORF dictionaries passing length compatibility
    :rtype: list[dict]
    """

    filtered = []

    uniprot_len = len(uniprot_protein_seq)

    # iterate over each orf dictionary to get amino acid length and check if it is within ±20% of reference uniprot length - if yes, then add into filtered
    for x in orfs:
        orf_aa_len = x.get("length_aa")
        # if this ORF has missing/invalid length_aa, move to the next ORF
        if orf_aa_len is None:
            continue
        if length_compatibility(orf_aa_len, uniprot_len):
            filtered.append(x)

    return filtered

def make_protein_aligner_blosum62() -> Align.PairwiseAligner:
    
    """
    Create a global protein sequence aligner with the BLOSUM62 matrix.
    This function initializes a Biopython PairwiseAligner for global protein alignment with the BLOSUM62 substitution matrix and
    standard gap penalties, and returns a configured aligner for pairwise alignment.
    :return: Biopython PairwiseAligner
    :rtype: Align.PairwiseAligner
    """

    # perform global sequence alignment
    aligner = Align.PairwiseAligner()
    aligner.mode = "global"

    # use BLOSUM62 (standard for protein sequence alignment)
    aligner.substitution_matrix = substitution_matrices.load("BLOSUM62")

    # gap penalties (standard for gap penalties)
    aligner.open_gap_score = -10
    aligner.extend_gap_score = -0.5
    return aligner

def seq_alignment_score(aligner: Align.PairwiseAligner, uniprot_aa: str, orf_aa: str) -> float:
    
    """
    Compute a global protein alignment score between UniProt amino acid sequence and candidate ORF amino acid sequence.
    This function uses a pre-configured Biopython PairwiseAligner to calculate a global alignment score between the UniProt
    reference protein and candidate ORF amino acid sequence.
    :param aligner: Pre-configured Biopython PairwiseAligner
    :type aligner: Align.PairwiseAligner
    :param uniprot_aa: UniProt reference protein amino acid sequence
    :type uniprot_aa: str
    :param orf_aa: Candidate ORF translated amino acid sequence
    :type orf_aa: str
    :return: Global alignment score
    :rtype: float
    """

    # return alignement score between UniProt amino acid sequence and candidate ORF amino acid sequence
    return aligner.score(uniprot_aa, orf_aa)


def select_best_orf_by_alignment(orfs: list[dict], uniprot_protein_seq: str, score_fraction_threshold=0.85) -> dict:

    """
    Identify the candidate ORF sequences that best matches a UniProt reference protein sequence.
    This function selects the best matching ORF from a list of candidate ORFs using a multi-step pipeline:
    1. Handle missing UniProt sequence or list of candidate ORFs
    2. Pre-filter ORFs by amino acid length (±20% of UniProt length)
    3. Compute global protein alignment scores using BLOSUM62 with gap penalties
    4. Select the ORF with the highest alignment score
    5. Determine whether the best ORF match passes a similarity threshold with UniProt reference protein sequence
    The final result is returned as a structured dictionary.
    :param orfs: List of candidate ORF dictionaries
    :type orfs: list[dict]
    :param uniprot_protein_seq: UniProt reference protein amino acid sequence
    :type uniprot_protein_seq: str
    :param score_fraction_threshold: Minimum threshold of the fraction for self-alignment score required to accept a match (default: 0.85)
    :type score_fraction_threshold: float
    :return: Dictionary with match status, ORF metadata, and alignment statistics
    :rtype: dict
    """

    # errors and checks for UniProt sequence and candidate ORF list
    if not uniprot_protein_seq or not isinstance(uniprot_protein_seq, str):
        return {
            "status": "no_match_found",
            "message": "UniProt protein sequence is missing or invalid."
        }
    if not orfs:
        return {
            "status": "no_orfs_found",
            "message": "ORFs are missing or invalid."
        }

    # filter by length for best ORF candidates - more efficient and gets rid of sequences that are too short or too long comapred to UniProt sequence length
    filtered_orf = filter_orf_length(orfs, uniprot_protein_seq)

    # if the length tolerance was too strict, fall back to scoring all candidate ORFs (unfiltered) in next step
    if not filtered_orf:
        filtered_orf = orfs

    # align and compare alignement scores
    # aligner is made once - prevents having to create aligner with BLOSUM62 over and over again through each iteration if inside for loop
    aligner = make_protein_aligner_blosum62()
    best_orf = None
    best_score = float("-inf")

    # iterate over each orf in the filtered list of ORFs to get alignment scores and to find the ORF with the best score
    for x in filtered_orf:
        orf_aa = x.get("orf_aa")
        if not orf_aa:
            continue
        score = seq_alignment_score(aligner, uniprot_protein_seq, orf_aa)
        if score > best_score:
            # best score and best_orf are the new highest alignment score and responsible ORF
            best_score = score
            best_orf = x

    if best_orf is None:
        return {
            "status": "no_match_found",
            "message": "ORFs were found, but could not be aligned."
        }

    # calculate score_fraction to see if the best ORF alignemnt score is a good enough match to the UniProt amino acid sequence
    # find the best maximum alignment score
    maximum_score = seq_alignment_score(aligner, uniprot_protein_seq, uniprot_protein_seq)

    if maximum_score <= 0:
        return {
            "status": "no_match_found",
            "message": "Unable to compute a valid maximum UniProt self-alignment score."
        }
    
    # calculate score fraction by normalising best ORF alignment score to maximum alignment score to determine if sequence similarity reached threshold
    score_fraction = best_score / maximum_score

    # copy ORF metadata and update with alignment statistics
    result = best_orf.copy()
    result.update({
        "alignment_score": best_score,
        "score_fraction": score_fraction,
        "filtered_orf_count": len(filtered_orf),
        "total_orf_count": len(orfs)
    })

    # outputs depending if score_fraction is equal to or above 0.85 threshold (match) or not (no match)
    if score_fraction >= score_fraction_threshold: # default set in function parameter is 0.85
        # print results, with updated status "match_found"
        result["status"] = "match_found"
    else:
        # print results, with updated status "no_match_found"
        result["status"] = "no_match_found"
        result["message"] = "The best ORF sequence does not meet the similarity threshold to the UniProt sequence."

    return result
