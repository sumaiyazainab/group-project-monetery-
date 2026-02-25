from Bio import SeqIO
from Bio.Seq import Seq

def validate_dna_sequence(sequence: str) -> None:

    """
    DNA sequence validation.
    This function validates DNA sequences and raises a ValueError if the sequence contains characters other than A, C, G, T (or N).
    :param sequence: DNA sequence
    :type sequence: str
    """

    seq = sequence.upper()
    nucleotides = {"A", "T", "C", "G", "N"}
    # make set instead of list - this gets rid of any dupliactes when printing, and sorts in ascending order
    inv_char = sorted({x for x in seq if x not in nucleotides})
#    for x in sequence: (instead of list comprehensions, this is another way - original for loop code I wrote)
#        if x not in nucleotides:
#            inv_char.append(x)
    if inv_char: #inv_char = [] is false, inv_char = [..] is true, therefore if inv_char contains elements, raise valueerror
        raise ValueError(f"Invalid characters found in DNA sequence: {inv_char}. Only A, C, G, T (and N) are allowed.")
    
def fasta_parsing(fasta_file) -> dict:

    """
    Parse plasmid FASTA file.
    This function parses a plasmid FASTA file, validates for appropriate DNA sequence, and extracts the header, nucleotide sequence, and sequence length.
    ValueError is raised if DNA sequence contains invalid characters that are not A, C, G, T (and N).
    :param fasta_file: FASTA file including the header and DNA sequence of one plasmid
    :type fasta_file: str
    :return: Dictionary of the plasmid header, DNA sequence, and length_bp.
    :rtype: dict
    """
    # SeqIO.parse accepts either a filepath or an open handle; made into a list to check for single-record plasmid sequence in file
    list_records = list(SeqIO.parse(fasta_file, "fasta"))
    # raise value error if there are more than one plasmid sequence in the uploaded fasta file
    if len(list_records) > 1:
        raise ValueError("There is more than one plasmid sequence in the FASTA file. Please check and try again.")
    # raise value error if there are no plasmid sequences inside the uploaded fasta file
    if len(list_records) == 0:
        raise ValueError("No plasmid sequence was found in the FASTA file. Please check and try again.")
    seq_record = list_records[0]
    sequence = str(seq_record.seq).upper()
    # validate DNA sequence before proceeding, so FLASK routes can surface a clean value error now before downstream pipeline
    validate_dna_sequence(sequence)
    plasmid_dict = {
        "header" :  seq_record.id,
        "sequence" : sequence,
        "length_bp" : len(sequence)
    }
    return plasmid_dict


def translate_with_alt_start(sequence: str, start_codons: set[str] | None = None) -> str:
   
    """
    This function translates a string of DNA nucleotide sequence into amino acid sequence, stopping before the stop codon, and forces
    the first amino acid to M if the first codon is in start_codons (e.g. ATG, GTG, TTG).
    :param sequence: DNA nucleotide sequence
    :type sequence: str
    :param start_codons: Set of codons treated as valid translation initiators. If None, defaults to common bacterial start codons {ATG, GTG, TTG}.
    :type start_codons: set[str] or None
    :return: protein amino acid sequence
    :rtype: str
    """

    seq = sequence.upper()
    # if user did not provide start_codons selection, use default start_codons
    if start_codons is None:
        start_codons = {"ATG", "GTG", "TTG"}
    # using biopython translate module, translate sequence up until the stop codon - do not include *
    aa_seq = str(Seq(seq).translate(to_stop=True))
    first_codon = seq[:3]
    # if the first_codon in aa_seq is in the selected start_codons, force the first codon into M for UniProt reference sequence alignment
    if first_codon in start_codons and aa_seq:
        aa_seq = "M" + aa_seq[1:]

    return aa_seq



def candidate_orf(plasmid_dict: dict, start_codons: set[str] | None = None) -> list[dict]:

    """
    Find and translate candidate ORF from plasmid FASTA file.
    This function identifies candidate ORFs in all six reading frames from a validated, circular, plasmid FASTA file
    in preparation for UniProt protein sequence matching. Candidate ORFs are translated, and respective 
    metadata are extracted as a list of dictionaries. Final translasted ORF amino acid sequence excludes stop codon,
    for matching UniProt protein sequences (which excludes terminal *). Any protein sequence shorter than 50 amino acid
    is excluded from the final list.
    :param plasmid_dict: FASTA file including the header and DNA sequence of one plasmid
    :type plasmid_dict: dict
    :param start_codons: Set of codons treated as valid translation initiators. If None, defaults to common bacterial start codons {ATG, GTG, TTG}.
    :type start_codons: set[str] or None
    :return: List of dictionaries of candidate ORFs and metadata including:
    - strand (str): '+' or '-' coding strand
    - frame (str): Reading frame (f1/f2/f3/r1/r2/r3)
    - start_nt (int): 0-based index of the first nucleotide of the ORF
    - end_nt (int): End coordinate of the ORF on the plasmid. For the forward strand, this is an end-exclusive index on the unrolled circular sequence.
    For the reverse strand, this is a mapped coordinate within the original plasmid range. The wraps_origin flag indicates whether the ORF crosses the plasmid origin.
    - wraps_origin (bool): True if ORF crosses the origin (wraparound) when represented on the circular plasmid
    - length_nt (int): Length of ORF - nucleotide sequence, including stop codon
    - orf_nt (str): ORF nucleotide sequence as read on the coding strand (includes stop codon)
    - length_aa (int): Length of translated ORF in amino acids (excluding stop codon)
    - orf_aa (str): Translated ORF amino acid sequence (excluding terminal '*')
    :rtype: list[dict]
    """

    # extract validated plasmid sequence (assumed uppercase DNA)
    sequence = plasmid_dict["sequence"]
    # duplicate sequence to accomodate for circularity and allow ORFs that cross the plasmid origin to be detected
    circular_seq = sequence + sequence
    L = len(sequence)
    # requirements for a candidate ORF
    # 50 aa is a common guidline in prokaryotic ORF screening; short peptides lack reliable UniProt matches
    min_aa_len = 50
    # default: common bacterial initiators - takes into account 2 more alternative start codons
    # this increases sensitivity but also increases the number of candidate ORFs; downstream length/alignment filters reduces noise
    if start_codons is None:
        start_codons = {"ATG", "GTG", "TTG"}
    stop_codons = {"TAA", "TAG", "TGA"}

    orfs = []

    # forward strand: frames f1, f2, f3
    for i in range(3):
        # offset sequence to represent the current reading frame (f1, f2, f3)
        frame_seq = circular_seq[i:]
        # scan the linear circular sequence in-frame (in 3s) in the reading frame
        for x in range(0,len(frame_seq)-2,3):
            codon = frame_seq[x:x+3]
            # identify for a start codon (ATG)
            if codon in start_codons:
                # reset ORF nucleotide sequence for each new start codon
                orf_nt = ""
                # extend ORF nucleotide sequence until an in-frame stop codon is reached
                for y in range(x, len(frame_seq)-2,3):
                    codon_2 = frame_seq[y:y+3]
                    orf_nt += codon_2
                    if codon_2 in stop_codons:
                        # compute and store the coordinates relative to the original plasmid
                        start_coord = i + x
                        stop_coord = i + y + 3
                        length_nt = len(orf_nt)
                        # the ORF must not exceed one full plasmid length
                        if start_coord < L and length_nt <= L:
                            # ORF wraps the origin if it extends past the original length
                            wraps_origin = stop_coord > L
                            # translate ORF nucleotide sequence to amino acid sequence and remove terminal stop codon
                            # in preparation of comparison with UniProt sequence (starts with M, does not contain *)
                            orf_aa = translate_with_alt_start(orf_nt, start_codons = start_codons)
                            # apply a minimum amino acid length filter of 50 to reduce search space and false positives from short random ORFs
                            if len(orf_aa) >= min_aa_len:
                                orfs.append({
                                    "strand" : "+",
                                    "frame" : f"f{i+1}",
                                    "start_nt" : start_coord,
                                    "end_nt" : stop_coord,
                                    "wraps_origin": wraps_origin,
                                    "length_nt" : length_nt,
                                    "orf_nt" : orf_nt,
                                    "length_aa" : len(orf_aa),
                                    "orf_aa" : orf_aa
                                })
                            # stop scanning this ORF after the first stop codon
                            break


    # reverse strand: frames r1, r2, r3
    circular_seq_list = list(circular_seq)
    # reverse the linear circular sequence and obtain its reverse complement
    circular_seq_list.reverse()
    seq_dict={'A':'T','T':'A','G':'C','C':'G','N':'N'}
    circular_seq_2 = "".join(seq_dict[x] for x in circular_seq_list)

    for i in range(3):
        # offset sequence to represent the current reading frame (r1, r2, r3)
        frame_seq = circular_seq_2[i:]
        # scan the linear circular sequence in-frame (in 3s) in the reading frame
        for x in range(0,len(frame_seq)-2,3):
            codon = frame_seq[x:x+3]
            # identify for a start codon (ATG) on the reverse strand
            if codon in start_codons:
                # reset ORF nucleotide sequence for each new start codon
                orf_nt = ""
                # extend ORF nucleotide sequence until an in-frame stop codon is reached
                for y in range(x, len(frame_seq)-2,3): # build sequence in 3s from ATG using coord
                    codon_2 = frame_seq[y:y+3]
                    orf_nt += codon_2
                    if codon_2 in stop_codons:
                        # same requirements as the forward strand scan:
                        # - the ORF start must occur within the first plasmid length (start_rev_coord < L)
                        # - the ORF must not be longer than one full plasmid (length_nt <= L)
                        start_rev_coord = i + x
                        stop_rev_coord = i + y + 3
                        length_nt = len(orf_nt)
                        if start_rev_coord < L and length_nt <= L:
                            # the modulo wraps the reverse-complement positions back into a single-plasmid range (0, L) instead of (0, 2L)
                            start_pos = start_rev_coord % L 
                            stop_pos = stop_rev_coord % L
                            # intuition: position p on the reverse-complement corresponds to (L - p) on the original plasmid
                            start_coord = (L - stop_pos) % L
                            stop_coord = (L - start_pos) % L
                            # reverse ORFs wrap origin if the stop position precedes the start position
                            wraps_origin = stop_coord <= start_coord
                            orf_aa = translate_with_alt_start(orf_nt, start_codons = start_codons)
                            if len(orf_aa) >= min_aa_len:
                                orfs.append({
                                    "strand" : "-",
                                    "frame" : f"r{i+1}",
                                    "start_nt" : start_coord,
                                    "end_nt" : stop_coord,
                                    "wraps_origin": wraps_origin,
                                    "length_nt" : length_nt,
                                    "orf_nt" : orf_nt,
                                    "length_aa" : len(orf_aa),
                                    "orf_aa" : orf_aa
                                })
                            # stop scanning this ORF after the first stop codon
                            break

    return orfs
            
