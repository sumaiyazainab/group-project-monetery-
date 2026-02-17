```python
# imports
```


```python
import pandas as pd
import numpy as np
from Bio import SeqIO
from Bio.Align import PairwiseAligner
```


```python
aligner = PairwiseAligner()
aligner.mode = "global"
aligner.match_score = 1
aligner.mismatch_score = 0
aligner.open_gap_score = -1
aligner.extend_gap_score = -0.5
```


```python
# load wt reference sequence
```


```python
def load_wt_sequence(fasta_path):
    try:
        record = SeqIO.read(fasta_path, "fasta")
        return str(record.seq)
    except Exception as e:
        raise ValueError(f"Could not load WT FASTA: {e}")
```


```python
WT_FASTA_PATH = "pET-28a_BSU_DNA_Pol_I_WT.fa"
wt_sequence = load_wt_sequence(WT_FASTA_PATH)
```


```python
# sequence validation 
```


```python
def validate_sequence(seq):
    if not isinstance(seq, str):
        raise TypeError("Sequence is not a string")
    if len(seq) == 0:
        raise ValueError("Empty sequence")
    if not set(seq.upper()).issubset({"A","C","G","T"}):
        raise ValueError(f"Invalid DNA characters: {seq}")
```


```python
# alignment reconstruction
```


```python
def align_sequences(seq1, seq2):
    alignment = aligner.align(seq1, seq2)[0]

    aligned1 = []
    aligned2 = []

    p1 = 0
    p2 = 0

    for (s1,e1), (s2,e2) in zip(alignment.aligned[0], alignment.aligned[1]):
        while p1 < s1:
            aligned1.append(seq1[p1])
            aligned2.append("-")
            p1 += 1
        while p2 < s2:
            aligned1.append("-")
            aligned2.append(seq2[p2])
            p2 += 1
        for i in range(e1 - s1):
            aligned1.append(seq1[s1+i])
            aligned2.append(seq2[s2+i])
        p1 = e1
        p2 = e2

    while p1 < len(seq1):
        aligned1.append(seq1[p1])
        aligned2.append("-")
        p1 += 1
    while p2 < len(seq2):
        aligned1.append("-")
        aligned2.append(seq2[p2])
        p2 += 1

    return "".join(aligned1), "".join(aligned2)
```


```python
# wt detection and mutation distance
```


```python
def mutation_distance(seq, wt_seq):
    a_wt, a_seq = align_sequences(wt_seq, seq)
    return sum(x != y for x, y in zip(a_wt, a_seq))


def is_sequence_wt(seq, wt_seq):
    return mutation_distance(seq, wt_seq) == 0
```


```python
# phase 6 main function
```


```python
def process_sequences_phase6(df, sequence_col, wt_seq):
    if sequence_col not in df.columns:
        raise KeyError(f"Sequence column '{sequence_col}' not found")

    df = df.copy()
    wt_flags = []
    mut_dists = []

    for i, seq in enumerate(df[sequence_col]):
        try:
            validate_sequence(seq)
            dist = mutation_distance(seq, wt_seq)
            mut_dists.append(dist)
            wt_flags.append(dist == 0)
        except Exception as e:
            raise RuntimeError(f"Row {i} failed: {e}")

    df["mutation_distance"] = mut_dists
    df["is_wt_sequence"] = wt_flags
    return df
```


```python
# activity score 
```


```python
def compute_activity_score(df, protein_col, dna_col):
    if protein_col not in df.columns:
        raise KeyError(f"Missing column {protein_col}")
    if dna_col not in df.columns:
        raise KeyError(f"Missing column {dna_col}")

    wt = df[df["is_wt_sequence"] == True]
    if len(wt) == 0:
        raise ValueError("No WT rows detected")

    protein_wt = wt[protein_col].mean()
    dna_wt = wt[dna_col].mean()

    if protein_wt == 0 or dna_wt == 0:
        raise ZeroDivisionError("WT mean is zero")

    df = df.copy()
    df["activity_score"] = (df[protein_col]/df[dna_col]) / (protein_wt/dna_wt)
    return df
```


```python
# global wt normalisation (trend plots)
```


```python
def add_global_wt_normalised_score(df, protein_col):
    wt = df[df["is_wt_sequence"] == True]
    global_wt_mean = wt[protein_col].mean()

    df = df.copy()
    df["global_wt_normalised"] = df[protein_col] / global_wt_mean
    return df
```


```python
# example
```


```python
df = pd.read_csv("DE_BSU_Pol_Batch_1.tsv", sep="\t")
```


```python
df = process_sequences_phase6(
    df,
    sequence_col="Assembled_DNA_Sequence",
    wt_seq=wt_sequence
)
```


```python
df = compute_activity_score(
    df,
    protein_col="Protein_Quantification_pg",
    dna_col="DNA_Quantification_fg"
)
```


```python
df = add_global_wt_normalised_score(
    df,
    protein_col="Protein_Quantification_pg"
)
```


```python
df[["mutation_distance","is_wt_sequence","activity_score","global_wt_normalised"]].head()
```
