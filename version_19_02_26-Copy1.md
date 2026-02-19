```python
# sequence processing (phase 6)
```


```python
# imports
```


```python
import pandas as pd
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Align import PairwiseAligner
```


```python
# load data (tsv or json)
```


```python
def load_dataset(path):
    if path.endswith(".tsv"):
        return pd.read_csv(path, sep="\t")
    elif path.endswith(".json"):
        return pd.read_json(path)
    else:
        raise ValueError("Unsupported file type")
```


```python
# load optional reference fasta 
```


```python
def load_reference_fasta(path=None):
    if path is None:
        return None
    record = SeqIO.read(path, "fasta")
    return str(record.seq)
```


```python
# identify wt rows from tsv 
```


```python
def detect_wt_rows(df):
    wt_mask = (
        (df["Directed_Evolution_Generation"] == 0) &
        (df["Parent_Plasmid_Variant"] == -1)
    )
    return df[wt_mask]
```


```python
# dna to rna to protein
```


```python
def safe_translate_dna(dna):
    if pd.isna(dna):
        return None
    
    dna = dna.upper().replace(" ", "").replace("\n", "")
    
    # trim to multiple of 3
    trim_len = len(dna) - (len(dna) % 3)
    dna = dna[:trim_len]
    
    if len(dna) == 0:
        return None
    
    try:
        return str(Seq(dna).translate(to_stop=False))
    except:
        return None
```


```python
def translate_variants(df):
    df = df.copy()
    df["protein_seq"] = df["Assembled_DNA_Sequence"].apply(safe_translate_dna)
    return df
```


```python
# choose reference protein
```


```python
def choose_reference_protein(df, ref_dna):
    if ref_dna:
        ref_protein = safe_translate_dna(ref_dna)
    else:
        wt_rows = df[df["Directed_Evolution_Generation"] == 0]
        ref_protein = wt_rows["protein_seq"].iloc[0]
    
    if ref_protein is None or len(ref_protein) == 0:
        raise ValueError("Reference protein sequence is empty.")
    
    return ref_protein
```


```python
# align and detect mutations
```


```python
from Bio.Align import PairwiseAligner

aligner = PairwiseAligner()
aligner.mode = "global"

def detect_mutations(var_seq, wt_seq):
    if not var_seq or not wt_seq:
        return []

    alignment = aligner.align(wt_seq, var_seq)[0]

    wt_aligned = alignment.target
    var_aligned = alignment.query

    muts = []
    pos = 0

    for w_char, v_char in zip(wt_aligned, var_aligned):
        if w_char != "-" and v_char != "-":
            pos += 1
            if w_char != v_char:
                muts.append(f"{w_char}{pos}{v_char}")

    return muts
```


```python
# synonymous vs non-synonymous
```


```python
def classify_mutations(dna_seq, wt_dna):
    muts = []
    for i in range(0, min(len(dna_seq), len(wt_dna)), 3):
        codon_var = dna_seq[i:i+3]
        codon_wt = wt_dna[i:i+3]
        if codon_var != codon_wt:
            aa_var = str(Seq(codon_var).translate())
            aa_wt = str(Seq(codon_wt).translate())
            if aa_var == aa_wt:
                muts.append("synonymous")
            else:
                muts.append("nonsynonymous")
    return muts
```


```python
# compute (multi-mode) baseline
```


```python
def get_baseline(df, mode="WT"):
    if mode == "WT":
        base = detect_wt_rows(df)
    elif mode == "generation0":
        base = df[df["Directed_Evolution_Generation"] == 0]
    elif mode == "min_mut":
        base = df[df["mutation_count"] == df["mutation_count"].min()]
    else:
        raise ValueError("Invalid baseline mode")

    if len(base) == 0:
        raise ValueError("No baseline rows detected")

    return (base["Protein_Quantification_pg"] /
            base["DNA_Quantification_fg"]).mean()
```


```python
# compute activity score
```


```python
def compute_activity_score(df, baseline_mode="WT"):
    df["raw_activity"] = (
        df["Protein_Quantification_pg"] /
        df["DNA_Quantification_fg"]
    )

    baseline = get_baseline(df, baseline_mode)

    df["Activity_Score"] = df["raw_activity"] / baseline
    return df
```


```python
# final pipeline function 
```


```python
def run_pipeline(data_path, fasta_path=None, baseline_mode="WT"):
    df = load_dataset(data_path)
    ref_dna = load_reference_fasta(fasta_path)
    df = translate_variants(df)
    ref_protein = choose_reference_protein(df, ref_dna)

    df = df[df["protein_seq"].notna()]
    df = df[df["protein_seq"].str.len() > 0]

    df["mutations"] = df["protein_seq"].apply(
        lambda s: detect_mutations(s, ref_protein)
    )
    df["mutation_count"] = df["mutations"].apply(len)

    df = compute_activity_score(df, baseline_mode)
    return df
```


```python
#
```


```python
#
```


```python
#
```


```python
# testing the pipeline
```


```python
df = run_pipeline(
    data_path="DE_BSU_Pol_Batch_1.tsv",
    fasta_path="pET-28a_BSU_DNA_Pol_I_WT.fa",
    baseline_mode="WT"
)

df.head()
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
      <th>mutation_count</th>
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
      <td>LVAFDAGKTTFRHGTFKEYKGGRQKTPPELSEQMPFIRELLDAYQI...</td>
      <td>[W1L, R2V, M3A, G4F, R5D, L7G, *8K, R9T, R10T,...</td>
      <td>2468</td>
      <td>0.101687</td>
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
      <td>LT*ASIFVMLVRGAEPMEKRQQRGLFTVPGLLLAFCSHVLSCVIP*...</td>
      <td>[W1L, R2T, M3*, G4A, R5S, A6I, L7F, *8V, R9M, ...</td>
      <td>2458</td>
      <td>0.082019</td>
      <td>0.806584</td>
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
      <td>VRRLSPSSSKSGAPFRVPI*CFTAPRPQKT*LG*WFT*WAIALIDG...</td>
      <td>[W1V, M3R, G4L, R5S, A6P, L7S, *8S, R9S, R10K,...</td>
      <td>2455</td>
      <td>0.076568</td>
      <td>0.752980</td>
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
      <td>DKL*PSPGAACVRGFHRHHRNARGSCGKAHQRGREAIHRCLPVHPR...</td>
      <td>[W1D, R2K, M3L, G4*, R5P, A6S, L7P, *8G, R9A, ...</td>
      <td>2450</td>
      <td>0.079278</td>
      <td>0.779623</td>
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
      <td>HRSRSRCLMSELTYINCVALTARFPVGKPVVPAALMNRPTRGERRF...</td>
      <td>[W1H, M3S, G4R, R5S, A6R, L7C, *8L, R9M, R10S,...</td>
      <td>2484</td>
      <td>0.076755</td>
      <td>0.754811</td>
    </tr>
  </tbody>
</table>
</div>




```python
df["protein_seq"].str.len().min()
```




    2607




```python
df["mutation_count"].describe()
```




    count     301.000000
    mean     2462.132890
    std        12.723559
    min      2418.000000
    25%      2454.000000
    50%      2461.000000
    75%      2471.000000
    max      2494.000000
    Name: mutation_count, dtype: float64


