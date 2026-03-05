IMPORTS
* pandas Ð 2.3.3
* logging Ð 0.5.1.2
* Bio.Seq import Seq
* Bio.Align import PairwiseAligner

PANDAS

Pandas is a data analysis library built around a table-like structure (rows and columns, like Excel/SQL) Ð the DataFrame. In the analysis pipeline of this project, pandas handles the management of the input data (load -> clean -> analyse -> export):
- The userÕs input of the DNA variant table (metadata) in either a .tsv or .json format; it loads the TSV or JSON experimental datasets into DataFrames. 
- For quality control: check that the required columns exist and remove rows with missing biological measurements. 
- To iterate over variants: loop over each plasmid variant, and analyse its DNA sequence
- Joining biological results back to the dataset Ð appending mutation analysis results (protein sequence, mutation counts, etc.) to the original dataset.
- To establish the userÕs output and export processed results
- The activity scoring math used to compute biochemical activity and normalise it to wild-type baselines  


LOGGING

Logging is PythonÕs built-in event recording system that writes timestamped messages to a file (rather than printing them to the screen). It records which variant is being processed, when an ORF cannot be found and where the pipeline may fail or skip data. This allows us to maintain a permanent record of whatÕs being done and how long it takes, rather than having silent failures. This audit trail is especially important for reproducibility and error diagnosis in scientific pipelines like this one. Having a file like this makes debugging long-running biological analyses easier Ð especially when managing a system that processes large, complex datasets, where we need to track how long parts of the system take to complete and the source of any failures. Here, users receive a file (called Òphase6.logÓ) once variant analysis is complete.


BIO.SEQ IMPORT SEQ

Bio.Seq comes from Biopython. Seq represents a biological sequence (DNA or protein) with built-in biological operations. In this case, Seq handles DNA-to-protein translation, ensuring biologically correct translation. (It also works together with aligned sequences after Bio.Align). It performs codon lookup, genetic-code translation, and handles stop codons, instead of manually mapping codons. Seq also analyses codon-level mutations and is used to distinguish synonymous from nonsynonymous mutations. It also detects truncation, which relies on translation producing Ò*Ó for stop codons.


BIO.ALIGN IMPORT PAIRWISEALIGNER

Bio.Align provides sequence alignment algorithms, and PairwiseAligner optimally aligns two sequences based on a scoring scheme. Detailed under the imports, the scoring scheme details a global alignment model that will align full-length proteins end to end, reward matches (+2), penalise mismatches (-1) and also penalise gaps (insertions/deletions). The alignment scoring defines how the algorithm determines the best alignment between the WT protein and the variant protein. It converts biological ideas into numbers, then chooses the alignment with the highest total score. A global alignment aligns the entire length of both proteins from start to end, so the algorithm is forced to account for all amino acids and all indels, not just a local matching region, which is useful because we expect the variant protein to be a mutated version of the same gene, not just a small fragment. A match score of +2 gives that score if two aligned amino acids are identical, which rewards conservation and pushes the algorithm to line up identical residues Ð so the same amino acid means itÕs likely the same codon and therefore is a good alignment. A mismatch score of -1 is given if two aligned amino acids differ, lightly penalising substitutions, so an amino acid substitution indicates a plausible mutation. The gap penalties (indels) are -5. This means that a score of -5 is assigned when a gap starts, which strongly discourages new gaps unless necessary. A gap extend of -1 is the cost of a continued gap (each extra gap costs -1), so one long indel is more likely than many small ones. The total alignment score is computed by summing these contributions, and, in comparison with a bad alignment, the higher-scoring alignment is chosen. This enables the algorithm to identify the most biologically plausible correspondence between WT and variant proteins. The downstream logic assumes that gaps indicate indels and that mismatches indicate substitutions. If scoring were wrong, we might misclassify substitutions as indels, misalign codons or detect fake truncations, so the scoring scheme directly affects mutation counts, synonymous vs nonsynonymous classification and truncation detection. 
This supports WT and variant proteins differing in length (indels), so the pipeline now translates both sequences into proteins, aligns them properly, and uses that alignment to map proteins. In protein alignment, it is used to identify the optimal alignment between the WT and variant proteins. It builds gapped protein strings so it can explicitly see substitutions, insertions, deletions. It also converts protein alignments into codon alignments, enabling synonymous and nonsynonymous classification, indel detection, and truncation detection, thereby improving mutation detection.    


QUALITY CONTROL

ÔESSENTIAL_FIELDSÕ is the basis of quality control that defines the minimum set of columns that must be present for the pipeline to work. The fields in detailed in the dictionary are required because later steps depend on them. ÔAssembed_DNA_SequenceÕ is needed for ORF finding and translation, ÔControlÕ is needed to identify WT rows, ÔDNA-Quantification_fgÕ and ÔProtein_Quantification_pgÕ are needed for activity scoring. 

The Ôparse_file(file_path)Õ function handles the file input only. It detects the file type by extension (.tsv or .json), loads it into a pandas DataFrame and rejects unsupported formats. It separates I/O concerns from biological logic, so the rest of the pipeline wonÕt worry about where the data is coming from, as it always receives a DataFrame. This improves reproducibility, extensibility, and makes responsibility clearer.

The Ôrun_quality_control(df)Õ function carries out two checks: column validation and row completeness filtering. Column validation ensures the dataset has all the required columns. If any of the columns are missing, the pipeline stops immediately Ð this prevents silent failures later on (e.g. trying to access a column that doesnÕt exist). Row completeness filtering keeps only rows with no missing values in the required columns and separates out rows that were rejected. So, we end up with a ÔcleanÕ data frame that is safe to analyse and a ÔrejectedÕ data frame that failed QC. This protects downstream steps from missing DNA sequences, quantification values, or control labels, which would disrupt biological interpretation.

The Ôparse_and)_qc(file_path)Õ function provides a single clean entry point. It loads the file, runs QC and returns the cleaned data (and rejected rows). This means that the main pipeline doesnÕt have to worry about file formats, separately call QC or remember the correct order, meaning that every data that enters the biological pipeline has passed QC.

Overall, the functions under quality control are responsible for loading the dataset, enforcing a required column schema, removing incomplete records, and ensuring that only structurally valid and biologically usable data enters the variant analysis pipeline.




PHASE 6 Ð VARIANT ANALYSIS

Ôbest_approx_match_pos(seq, seed, max_mismatches = 20)Õ

The aim of this function is to find where a short ÔseedÕ sequence best matches inside a longer DNA sequence, even if there are some mismatches. This is used to locate the gene (ORF) inside a circular plasmid. It slides a window the same length as the seed along the sequence, counts mismatches at each position, and keeps track of the position with the fewest mismatches. If the mismatch count is less than or equal to the Ômax_mismatchesÕ, it returns (best_position, mismatch_count); otherwise, it returns (none, best_mismatch_count). Mutated genes will not perfectly match the reference gene, so instead of exact matching, this allows substitutions, mild divergence and robust ORF detection. 

Ôextract_cds_to_first_stop(circular_seq:str, start_pos:int, max_cds_nt:int = 20000) -> strÕ

This extracts a CDS from a circular plasmid, starting at a given position and stopping at the first in-frame stop codon. So, weÕre reading forward till translation should stop from a predicted start codon. This treats the plasmid as circular (so indexing wraps around the end of the string), reads the codons starting at Ôstart_posÕ, continues appending codons until a stop codon is encountered or Ômax_cds_ntÕ is reached (the safety limit) and returns the nucleotide string of the CDS. Because your gene is embedded inside a plasmid, plasmids are circular, and the gene may start near the end of the sequence. This function ensures the ORF is biologically valid (starts at a start codon and ends at a stop codon).

Ôfind_orf_by_seed_and_stop(variant_seq, best_orf, seed_len=50, max_mismatches=20)Õ

This finds the full ORF within a plasmid by locating an approximate match to a known reference ORF (Ôthe seedÕ) and extracting the CDS from that position to the stop codon, identifying that gene in a mutated plasmid. The function takes the plasmid sequence and the known reference ORF and uses the first Ôseed_lenÕ bases of the reference ORF as a search probe. It then calls Ôbest_approx_match_pos()Õ to find where that seed best matches within the plasmid (allowing mismatches). When the best match is found, that position is treated as the start of the gene and Ôextract_cds_to_first_stop()Õ is used to extract the ORF from that position to the first stop codon. Because exact matching fails when mutations exist, this function allows substitutions and divergence and still finds the gene Ð without this, any mutation in the first part of the gene could break ORF detection. 

Ôtranslate_with_alt_start_with_stop(sequence:str, start_codon:set[str] | None=None -> strÕ

This translates the DNA into protein whilst allowing alternative start codons and preserving internal stop codons, allowing translation that mimics a real ribosome, not a nave string converter. It defines acceptable start codons, translates codons in triplets using the genetic code, doesnÕt truncate at the first stop codon (internal * characters are preserved and used later detect truncations) and returns a full amino acid string. Whereas standard translation functions often stop at the first * and assume only ÔATGÕ is valid, biologically bacteria use alternative start codons, and premature stops are meaningful mutations Ð this function preserves that information.

Ôglobal_align_proteins(wt_prot;str, var_prot:str)

The purpose of this is to perform a global alignment between WT and variant protein sequences (it lines up the two proteins so that corresponding residues match). It uses a global alignment algorithm (Needleman-Wunsch), applies the scoring schema (discussed earlier) and produces two aligned strings. Indels destroy positional correspondence, so without alignment, codon 50 wouldnÕt equal codon 50 in the variant protein Ð this alignment restores the homologous positions between proteins.




Ôbuild_codon_aligned_cds(wt_aln:str, var_aln:str, wt_nt_orf:str, var_nt_orf:str)Õ

This converts the protein alignment back into a codon-aligned DNA alignment. So, if thereÕs a gap in the protein space, this inserts a triple gap in the DNA space. The function walks through the aligned protein strings, keeping two pointers (one for WT DNA codons and one for variant DNA codons). And for each aligned position: if the protein = AA, it appends to the next codon; if the protein is = Ô-Ô, it appends Ô---Ô. Mutation classification occurs at the codon level, not just the protein level, so this step ensures that codons stay in-frame, that gaps are biologically meaningful, and that synonymous vs nonsynonymous mutations can be detected.

Ôclassify_muatiosn_codon_aligned(wt_cds:str, var_cds:str)Õ

The purpose of this is to count and classify mutations using aligned codon sequences. This function goes through codons and each codon pair it: checks if both codons are equal (no mutation in that case), checks if theyÕre different in which case, it translates each codon and compares the amino acids, classifies the mutations (synonymous, nonsynonymous, insertion (WT codon is Ô---Ô), deletion (variant codon is Ô---Ô)) and also detects premature stop codons (truncation). This is where biological interpretation comes in, and we try to answer key questions. Has any change to the protein occurred? If any, then how many? Was it truncated? Were the residues indels? These are key contributors to the variant analysis pipeline's output.

Ôextract_mutation_annotaions(wt_aln, var_aln)Õ

The purpose of this is to generate mutation labels that we can read from two aligned protein sequences (e.g. ÔT599SÕ, ÔA600delÕ, etc). This walks through the aligned proteins position by position and tracks the true WT position (Ôwt_posÕ), ignoring gaps. Then it classifies the mutations and produces mutation labels, converting raw alignment strings into biologically interpretable mutation notation.

Ôvariant_analysis(data-path, best_orf, seed_len=80, max_mismatches=8)Õ

This is the main pipeline function which coordinates the entire mutation analysis. 

1. Load and QC data, which ensures that the required columns exist and that there are no missing essential values, thus preventing broken biological analysis
2. Identify WT plasmid Ð the control plasmid carries the WT gene 
3. Extract the WT ORF. This uses the seed-based approximate match algorithm to locate the gene Ð if this fails, the whole thing aborts 
4. Translate WT gene, which creates the reference protein for alignment
5. Loop over all the variants and analyse one mutant plasmid
6. Extract variant ORF. If this fails, the warning is logged, it outputs empty results and continues to the next variant Ð this prevents the pipeline from crashing
7. Translate the variant gene from DNA to protein
8. Align the WT and variant proteins. This inserts gaps, aligns indels and restores positional correspondence, which is essential for correct mutation calling
9. Generate mutation labels and produce mutation strings and affected positions
10. Build codon-aligned DNA sequences. Turn protein gaps to DNA triplet gaps, which preserves the reading frame for codon-level analysis
11. Classify mutations. This computes the total mutations, synonymous vs nonsynonymous, insertion count, deletion count, and truncation flag. This helps us see whether the protein changed, by how much, and whether it was truncated.
12. Store the results. So each variant yields: protein_seq, mutations, aa_positions, mutation_count, synonymous, nonsynonymous, insertions, deletions and truncating
13. Merge with the original dataset Ð attaches the mutation data as new columns
14. Save outputs as .csv and .json. These are the final products of this pipeline.


DEVELOPING THE ACTIVITY SCORE METRIC 

The experiment measures DNA output (i.e., how much plasmid is present) and protein output. The analysis treats DNA abundance as the selected signal Ð but itÕll be corrected for protein amount to avoid confounding by expression/loading.

* Raw selection signal = DNA produced per unit protein (the pipeline uses protein as a denominator to control for expression/loading differences)
* Metric = High DNA output per unit protein
* Normalise (Ôraw_activityÕ): 
o Efficiency = DNA_yield / Protein_yield
* Then baseline correct using WT: 
o Activity Score = (DNA_variant / Protein_variant) / (DNA_WT / Protein_WT)
* This means: 
o WT Å 1.0
o Better than WT > 1
o Worse than WT < 1
* WT becomes the reference Å 1
* Variants scale relative to it 

What is the Activity Score trying to measure?
We want to estimate the intrinsic catalytic efficiency of a polymerase variant (not the total DNA made, total fluorescence, or total protein expressed); how much DNA the enzyme produces per unit enzyme.
Efficiency Å Product formed/Enzyme present
       Product = dsDNA (fluorescence)
       Enzyme = polymerase (A280 absorbance)
Therefore, the Activity Score is a proxy for catalytic efficiency rather than raw productivity.
Raw fluorescence is scientifically flawed because if selection is done on the fact that fluorescence is proportional to the DNA amount, it confounds two variables: 
	DNA yield = Ä(enzyme efficiency, enzyme concentration)
This means that a variant that is bad but highly expressed may look good, and a variant that is excellent but poorly expressed may look bad Ð a confounding variable problem. Statistically, this means that the protein concentration is a nuisance variable, and the catalytic efficiency is the latent variable of interest. The Activity Score must remove this confounder.
A good Activity Score should achieve the following criteria:
1. Biology interpretability Ð it canÕt be just an arbitrary number; it needs to correspond to something biologically meaningful. This allows us to make comparisons like ÒVariant X is 3 times more efficient than the WTÓ. 
2. Normalisation across experiments and generations. The absolute fluorescence values depend on: 
* Reagent concentrations
* Instrument settings
* Batch effects
* Temperature
* Dye lots
Therefore, the score must be relative rather than absolute, which requires normalisation to WT controls.
3. Penalise trivial solutions Ð the metric must penalise variants that only succeed by overexpression and reward variants that produce more DNA with less protein. Otherwise, weÕd select expression artefacts instead of better enzymes.
4. Robustness to noise. The experimental readouts are noisy due to droplet variation, fluorescence noise, and absorbance noise. Therefore, a good metric should: average WT baselines, avoid division by very small numbers, and reduce sensitivity to outliers.
5. Comparability across generations Ð it should be valid to plot Ôactivity score vs generationÕ so generation 1 and generation 8 scores are comparable. This requires a consistent baseline and consistent scaling. 
The ratio model is our core idea:
Ei = Di / Pi 
	Di = DNA yield variant i 
	Pi = protein yield of variant i 
This estimates the DNA produced per unit protein, but it suffers from scale differences between experiments, unknown units and batch effects.
The WT-normalised ratio, however, is an ideal dimensionless relative efficiency score. 
First, we define the WT efficiency:
	EWT = <DWT  / PWT >
Then: 
	Activity Scorei = (Di  / Pi)/EWT
This way, we achieve the points set out in the beginning, where:
o WT Å 1.0
o Variant > 1 = Better than WT
o Variant < 1 = Worse than WT
  This gives us the interpretability, comparability, batch correction, and fairness that weÕre aiming for. 
Generation-specific vs Global baseline
Global:
EWT = mean of all WT controls
Pros: 
- Stable
- Comparable across generations
Cons:
- Assumes conditions constant

Per-generation baseline:
EWT,g
Pros:
- Controls for drift
- More precise
Cons:
- Complicates cross-generation comparison

Compromise:
- Compute both
- Use per-generation for ranking
- Global for plotting trends

The global baseline is the mean of all ÔControl == TrueÕ rows Ð rather just Ôgen0Õ, the first row in the data frame. This treats the WT not as a single construct, but as a distribution of baseline measurements across the experiment. So, weÕre looking at how each variant compares to the typical parental construct under experimental conditions. This means that the WT is re-measured in each generation, the conditions drift, noise accumulates, and controls exist specifically to estimate the baseline. Biologically, each generationÕs control remains WT relative to its parent, and globally, they approximate the ancestral functional level. This reduces noise, makes it more stable, more reproducible, and robust to a single bad normalisation. Our dataset is from directed evolution, not a single-timepoint experiment. Each generation has its own parent, its own control, and its own experimental context, so the WT is not biologically frozen; itÕs the best-known parent at that generation. This means that using only generation 0 ignores assay drift, using all control pools information, and using per-generation controls captures stepwise evolution. Therefore, using the mean of all control measurements as a global baseline provides a lower-variance estimator of WT activity rather than relying on a single ancestral measurement, while retaining a biologically interpretable reference scale for cross-generation comparison. 
The generation-specific baseline (local reference) is the variants compared only to their immediate parent. This is good for selection efficiency per round, mutational-step effects, and local fitness landscapes. WeÕre now modelling the mutant vs its evolutionary parent, which is more faithful to directed evolution theory.

LIMITATIONS
* Division instability
       If Pi  -> 0, then Di  / Pi -> °
Biologically, an extremely low protein count yields an unreliable signal, so we should have imposed a minimum protein threshold or added a pseudo-count (Di / Pi + ?).
* Distribution skew 
Since ratios tend to be heavy tailed. A possible correction could have been log(Di / Pi), and then normalise it to the WT, Scorei = log(Di / Pi) Ð log(EWT). This gives symmetry, better statistics and easier plotting. 
* Statistical framing 
Frame Activity Score as:
	Observed DNA = ?. Protein + ?
Then estimate:
	?i = Di  / Pi 
That way, the Activity Score Å slope of DNA vs protein, variant specific. This connects our metric to linear regression, enzyme kinetics and efficiency estimation. 

PHASE 7 Ð ACTIVITY SCORING
Ôidentify_wt_rows(df)Õ
The purpose is to identify the rows in the dataset that correspond to WT (baseline) controls. It selects the rows where ÔControl == TrueÕ so the WT rows are not inferred from the generation number, they are explicitly labelled in the data. If no such rows exist, it raises an error. In this experiment, each generation has its own ÔparentÕ construct, and WT is defined experimentally rather than mathematically. So, this function encodes the biological WT as whatever the experimenter marked as ÔControlÕ.
Ôcompute_raw_activity(df)Õ
This function computes the uncorrected biochemical activity metric for each row Ð based on the earlier definition, Ôraw_activity = DNA_Quantification_fg / Protein_Quantification_pgÕ. This estimates the amount of DNA output per unit of protein and corrects for loading and expression differences, as well as total protein abundance. So the raw DNA alone is not to be trusted and needs to be normalised by protein. This returns a copy of the DataFrame with a new column Ôdf[Ôraw_activityÕ]Õ.
Ôcompute_global_wt_baseline(df)Õ
This computes a single WT reference value using all WT rows. It calls Ôidentify_wt_rows(df)Õ, takes their raw activity and computes the mean. This represents the average activity of all parental constructs across the experiment, assuming that the WT is stable and that a single global reference is valid. This is useful when weÕre comparing across generations and want a single fixed reference scale. 

Ôcompute_generation_wt_baselines(df)Õ
This function computes separate WT baselines for each generation. It filters WT rows (ÔControl == TrueÕ) and groups them by ÔDirected_Evolution_GenerationÕ and then computes the mean Ôraw_activityÕ for each group. 
This results in something like:
	Generation 0 ? baseline = 0.95  
       Generation 1 ? baseline = 1.02  
       Generation 2 ? baseline = 0.98
This accounts for drift between generations, batch effects, and changing parental backgrounds, so instead of a single WT baseline, we get a local WT baseline per generation. 
compute_activity_score(df, baseline_mode="global")
The purpose of the function is to compute the final normalised activity score for each variant. This is where normalisation happens. The user chooses a baseline strategy, ÔglobalÕ and ÔgenerationÕ, and it uses their respective functions (detailed above). It then computes the normalised score and adds a column, Ôdf[Ôactivity_scoreÕ] Ð and this is the number used for downstream interpretation. 
activity_scoring(df, baseline_mode="global")
This function orchestrates phase 7, tying all the steps above together in the appropriate order. It computes the raw activity, computes the normalised activity score, and returns the scored DataFrame (turning DNA, Protein into Ôraw_activityÕ + Ôactivity_scoreÕ) in one call.
