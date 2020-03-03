# Software develped or optimised for Papareddy et al. (2020)!

# Small RNA analysis!
Nextflow pipe line to analysed small RNA high throutput data sets made by compailing custom software used developed in Plotnikova et al. (2019).

General Description: Cutadapt (Martin, 2011) was used to trim adapter sequences from sRNA-seq reads and 18–30-base sequences that contained an adapter were retained. The trimmed sequences were aligned to the Arabidopsis thaliana TAIR10 genome (Lamesch et al., 2012) with STAR (Dobin et al., 2013) requiring no mismatches and allowing ≤100 multiple end-to-end alignments. Resulting SAM files were then processed with the readmapIO.py script to re-assign multimappers with a “rich-get-richer” algorithm as previously described (Schon et al., 2018). Output bedFiles were sorted, condensed and normalized for total genome-matching reads. The BEDtools map function (Quinlan and Hall, 2010) was then used to quantify sum of the normalised Reads per million mapping to TAIR10 annotated Transposable elements (TEs). Statistical analyses and associated figures were generated with the R statistical computing package (R Core Team, 2018).

Required softwares: You will need the following software: 1. Cutadapt (cutadapt/1.18-foss-2018b-Python-2.7.15) 2. STAR (star/2.7.1a-foss-2018b) 3.SAMTools (samtools/1.9-foss-2018b) 3.deepTools (deeptools/3.1.2-foss-2018b-python-2.7.15) 4.readmapIO.py (included, from Schon et al. (2018) Genome Research) 5. fasta_utils.py (included, from Schon et al. (2018) Genome Research) 6. BEDtools (bedtools/2.27.1-foss-2018b) 7. bed_collapse.py (included, this study) 8. FastQC (fastqc/0.11.8-java-1.8)

And the following files: 1. FASTA file of reference genome (not included but can be downloaded from ftp://ftp.ensemblgenomes.org/pub/plants/release- 44/fasta/arabidopsis_thaliana/dna/) 3. BEDfile of TAIR10 TEs thats are classified as class A/Euchromatic or class B/Heterochromatic(TEs.siRNA.clusterID.sorted)



This procedure was performed for the 25 sRNA-seq datasets generated as part of this study and all the analysed public;y availabile datasets described in sheet1 of dataset Table_S1_Mapping_stats.xlxs. 
