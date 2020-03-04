# Software develped or optimised for Papareddy et al. (2020)

## Small RNA analysis
Nextflow pipe line to analysed small RNA high throutput data sets; Made by compailing custom software pieces used in Plotnikova et al. (2019).

General Description: Cutadapt (Martin, 2011) was used to trim adapter sequences from sRNA-seq reads and 18–30-base sequences that contained an adapter were retained. The trimmed sequences were aligned to the Arabidopsis thaliana TAIR10 genome (Lamesch et al., 2012) with STAR (Dobin et al., 2013) requiring no mismatches and allowing ≤100 multiple end-to-end alignments. Resulting SAM files were then processed with the readmapIO.py script to re-assign multimappers with a “rich-get-richer” algorithm as previously described (Schon et al., 2018). Output bedFiles were sorted, condensed and normalized for total genome-matching reads. The BEDtools map function (Quinlan and Hall, 2010) was then used to quantify sum of the normalised Reads per million mapping to TAIR10 annotated Transposable elements (TEs). Statistical analyses and associated figures were generated with the R statistical computing package (R Core Team, 2018).

Required softwares: You will need the following software: 
1. Cutadapt (cutadapt/1.18-foss-2018b-Python-2.7.15) 
2. STAR (star/2.7.1a-foss-2018b) 
3. SAMTools (samtools/1.9-foss-2018b) 
4. deepTools (deeptools/3.1.2-foss-2018b-python-2.7.15) 
5. readmapIO.py (included, from Schon et al. (2018) Genome Research) 
6. fasta_utils.py (included, from Schon et al. (2018) Genome Research) 
7. BEDtools (bedtools/2.27.1-foss-2018b) 
8. bed_collapse.py (included, this study) 
9. FastQC (fastqc/0.11.8-java-1.8) and  
9. nextflow (nextflow/19.10.0)

And the following files: 1. FASTA file of reference genome (ftp://ftp.ensemblgenomes.org/pub/plants/release- 44/fasta/arabidopsis_thaliana/dna/) 2. BEDfile of TAIR10 TEs thats are classified as class A/Euchromatic or class B/Heterochromatic(Not included yet, will be released soon or on request).<br/>
<br/>
Dependency scripts and files for sRNA analyses as part of nextflow from our lab are in the folder sRNA_assets.

This procedure was performed for the 25 sRNA-seq datasets generated as part of this study and all the analysed publicly availabile datasets described in sheet1 of dataset Table_S1_Mapping_stats.xlxs. 

### Running nextflow pipeline
This pipline can be exicuted by running the command: <br/>
__nextflow run sRNA.nf__ <br/>
* Default setings are as described in the mauscripts but parameters can be changes according to the users' wish (see main.nf parameters section) ; *for example changing input file path, output file path and desired 3 prime adaptor sequences to be trimmed are as follows*<br/>
__--files__ 'PATH/TO/FASTQ/ <br/>
__--output__ 'PATH/TO/OUTPT/RESULTS/' <br/>
__--adapt__ 'ADAPTORSEQUENCE' <br/>

## MethylC-Seq analysis
General Description: Sequenced reads were trimmed using Trim Galore with default settings. In addition, the first six bases of each read were removed to exclude the random hexamer portion of the read used during the preamplification step. After quality filtering and adaptor trimming, bisulfite-converted reads were aligned against the TAIR10 genome (Lamesch et al., 2012) using Bismark (bismark --non_directional -q --score-min L,0,-0.4) (Krueger and Andrews, 2011). BAM files containing clonalonly deduplicated and uniquely mappedping reads were then used as input for the Methylpy software (https://bitbucket.org/schultzmattd/methylpy) to extract weighted methylation rate at each cytosine as previously described (Schultz et al., 2015). Bisulfite conversion rates were calculated using the unmethylated chloroplast genome or spiked-in unmethylated Lambda phage DNA (European Nucleotide Archive Accession Number J02459, Promega catalog number D1521) controls.
Differentially methylated regions (DMRs) were defined using Methylpy as described (Kawakatsu et al., 2017). Briefly, biological replicates were pooled and differentially methylated sites (DMSs) were identified by the root mean square tests with false discovery rates ≤ 0.01. Cytosine sites with ≥4 overlapping reads were retained for all samples except for preglobular in which DMSs with ≥3 overlapping reads were retained. Differentially methylated sites within 100-bp were collapsed into DMRs. CHH-DMRs were further filtered by discarding regions with < 4 DMSs and methylation difference < 20%. Using these parameters, DMRs were identified in all pairwise combinations across samples of interest and merged using the BEDtools merge function. DMRs were used to calculate the weighted CHH methylation rate on all analyzed tissue types.

__Nextflow pipeline combining Bismark (for MethylC-seq read Mapping) and Methylpy (for DMR calling and other useful miscellaneous) as done in this manuscript is under construction__<br/>
Meanwhile you look into following Git pages for details<br/>
<br/>
* https://github.com/FelixKrueger/TrimGalore
* https://github.com/FelixKrueger/Bismark
* https://github.com/yupenghe/methylpy
