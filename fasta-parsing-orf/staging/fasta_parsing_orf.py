from Bio import SeqIO

def validate_dna_sequence(sequence: str) -> None:

    """
    DNA sequence validation.
    This function validates DNA sequences and raises a ValueError if the sequence contains characters other than A, C, G, T (or N).
    :param sequence: DNA sequence
    :type sequence: str
    """

    nucleotides = {"A", "T", "C", "G", "N"}
    # this gets rid of any dupliactes when printing, and sorts in ascending order
    inv_char = sorted({x for x in sequence if x not in nucleotides})
    if inv_char:
        #if sequence contains any invalid characters, raise explicit valueerror
        raise ValueError(f"Invalid characters found in DNA sequence: {inv_char}. Only A, C, G, T (and N) are allowed.")
    

def fasta_parsing(fasta_file) -> dict:

    """
    Parse plasmid FASTA file.
    This function parses a plasmid FASTA file, validates for appropriate DNA sequence, and extracts the header, nucleotide sequence, sequence length and GC percentage.
    ValueError is raised if DNA sequence contains invalid characters that are not A, C, G, T (and N).
    :param fasta_file: FASTA file including the header and DNA sequence of one plasmid
    :type fasta_file: str
    :return: Dictionary of the plasmid header, DNA sequence, length_bp, and gc_percent.
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

def translate(sequence: str, start_codons: set[str] | None = None) -> str:
    
    """
    Translate DNA nucleotide sequence into protein amino acid sequence.
    This function translates a string of nucleotide sequence into amino acid sequence, including stop codon.
    :param sequence: DNA nucleotide sequence
    :type sequence: str
    :param start_codons: Set of codons treated as valid translation initiators. If None, defaults to common bacterial start codons {ATG, GTG, TTG}.
    :type start_codons: set[str] or None
    :return: protein amino acid sequence
    :rtype: str
    """
    codon_table = {
        'ATA':'I', 'ATC':'I', 'ATT':'I', 'ATG':'M',
        'ACA':'T', 'ACC':'T', 'ACG':'T', 'ACT':'T',
        'AAC':'N', 'AAT':'N', 'AAA':'K', 'AAG':'K',
        'AGC':'S', 'AGT':'S', 'AGA':'R', 'AGG':'R',
        'CTA':'L', 'CTC':'L', 'CTG':'L', 'CTT':'L',
        'CCA':'P', 'CCC':'P', 'CCG':'P', 'CCT':'P',
        'CAC':'H', 'CAT':'H', 'CAA':'Q', 'CAG':'Q',
        'CGA':'R', 'CGC':'R', 'CGG':'R', 'CGT':'R',
        'GTA':'V', 'GTC':'V', 'GTG':'V', 'GTT':'V',
        'GCA':'A', 'GCC':'A', 'GCG':'A', 'GCT':'A',
        'GAC':'D', 'GAT':'D', 'GAA':'E', 'GAG':'E',
        'GGA':'G', 'GGC':'G', 'GGG':'G', 'GGT':'G',
        'TCA':'S', 'TCC':'S', 'TCG':'S', 'TCT':'S',
        'TTC':'F', 'TTT':'F', 'TTA':'L', 'TTG':'L',
        'TAC':'Y', 'TAT':'Y', 'TAA':'*', 'TAG':'*',
        'TGC':'C', 'TGT':'C', 'TGA':'*', 'TGG':'W',
    }

    # if user doesn't specify the start_codon mode, then just use default
    if start_codons is None:
        start_codons = {"ATG", "GTG", "TTG"}

    aminoacid = ""
    
    # loops through the sequences for every 3 nucleotides, for each codon
    for x in range(0,len(sequence),3): 
        codon = sequence[x:x+3]
        if x == 0 and codon in start_codons:
            # X is used here for invalid codons/codons not found in dict
            aminoacid += "M"
        else:
            # if not a start codon, get the corresponding amino acid symbol
            # X is used here for invalid codons/codons not found in dict
            aminoacid += codon_table.get(codon,"X")

    return aminoacid




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
                            # in preparation of comparison with UniProt sequence (does not contain *)
                            orf_aa = translate(orf_nt, start_codons = start_codons).rstrip("*")
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
                            orf_aa = translate(orf_nt, start_codons = start_codons).rstrip("*")
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
            
            
