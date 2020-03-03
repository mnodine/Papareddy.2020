#!/bin/bash
#Pipeline to process sRNA-seq data, align to Col-0 genome and quantify number of reads per annotated miRNA loci
#Begin PBS directives
#PBS -P rnaseq_nod
#PBS -N sRNA.pipe
##PBS -J 1
#PBS -j oe
#PBS -o logs/sRNA.pipe/sRNA.pipe.log
#PBS -l walltime=10:00:00
#PBS -l select=1:ncpus=8:mem=32gb

#File paths to various directories
DIR=test										##Base name of directory
sRNA_ROOT=/lustre/scratch/users/michael.nodine/seq/sRNA/$DIR/sRNAseq/     		##Where sRNA-seq data is located
ALIGN_ROOT=/lustre/scratch/users/michael.nodine/seq/sRNA/$DIR/alignments/STAR_100/	##Where sRNA-seq alignments should be stored
BED_ROOT=/lustre/scratch/users/michael.nodine/seq/sRNA/$DIR/bedFiles/STAR_100/		##Where sRNA-seq alignments as bedFiles should be stored
SCRIPTS=/lustre/scratch/users/michael.nodine/seq/sRNA/$DIR/scripts/			##Where python scripts are located
GENOME=Col_0										##Name of genome aligning to
GENOME_DIR=/lustre/scratch/users/michael.nodine/compFiles/genomes/Ath/Col_0		##Where FASTA and indexed genome files for STAR aligments are located
SUFFIX=adaptertrimmed.fastq 								##File ending for FASTQ of adapter-trimmed sequences
MISMATCH=0										##Number of mismatches to allow during aligments
NAME=test										##Name of sample
FRAG_OVERLAP=0.8									##Proportion of read that must overlap annotated miRNA locus and vice-versa
ANNO_ROOT=/lustre/scratch/users/michael.nodine/seq/sRNA/$DIR/annotations/		##Where annotation files are located
ANNO=miRNA_mature									##Name annotation bedFile 

#1. Use Cutadapt to remove 3' adapters, and only keep reads that had 3' adapter and were at 18-30 bases long
##Sequence of 3' adapter used in NEBnext kit: AGATCGGAAGA

module load cutadapt/1.18-foss-2018b-Python-2.7.15
cutadapt -a AGATCGGAAGA -o $sRNA_ROOT/${NAME}_${SUFFIX} -m 18 -M 30 --trimmed-only $sRNA_ROOT/${NAME}.fastq

#2. Align to Col-0 TAIR10 genome with STAR allowing up to 100 multiple end-to-end alignments
module load rna-star/2.5.2a-foss-2016a
##First have to generate genome index if this hasn't been done already
STAR --runThreadN 4 --runMode genomeGenerate --genomeDir $GENOME_DIR --genomeFastaFiles $GENOME_DIR/${GENOME}.fa 

##Make directory for STAR output
mkdir -p $ALIGN_ROOT/${NAME}/

STAR --runThreadN 8 --genomeDir $GENOME_DIR --readFilesIn $sRNA_ROOT/${NAME}_${SUFFIX} --outFilterMultimapNmax 100 --outSAMtype SAM --outFilterMismatchNmax $MISMATCH --alignEndsType EndToEnd --outReadsUnmapped Fastx --alignIntronMax 1 --outFileNamePrefix $ALIGN_ROOT/${NAME}/ --outSAMprimaryFlag AllBestScore

#3. Re-assign multi-mappers using "rich-get-richer" algorithm
mkdir -p $BED_ROOT/${NAME}/
python $SCRIPTS/readmapIO.py -S $BED_ROOT/${NAME}/${NAME} -I $ALIGN_ROOT/${NAME}/Aligned.out.sam --bed -F $GENOME_DIR/${GENOME}.fa  --allow_naive --untemp_out ACGTN 

#4. Sort resulting BED files and collapse overlapping reads
module load BEDTools/2.27.1-GCCcore-6.4.0
bedtools sort -i $BED_ROOT/${NAME}/${NAME}.BODY.temp.bed > $BED_ROOT/${NAME}/${NAME}.sorted.bed
python $SCRIPTS/bed_collapse.py $BED_ROOT/${NAME}/${NAME}.sorted.bed > $BED_ROOT/${NAME}/${NAME}.bed

#5. Make normalized version of bedfile
##get total number of genome-mapping reads from BG file and use it to normalize BG hit-normalized reads
GMR=$(awk -v OFS='\t' '{sum+=$5}END{print sum}' $BED_ROOT/${NAME}/${NAME}.bed)
MGMR=$(awk '{print +$1 }' <<< $GMR)
MGMR=$(bc -l <<< $MGMR/1000000)
echo $MGMR

awk -v OFS='\t' -v MGMR="$MGMR" '{print $1,$2,$3,$4,$5/MGMR,$6,$7,$8}' $BED_ROOT/${NAME}/${NAME}.bed > $BED_ROOT/${NAME}/${NAME}.norm.bed

#6. Generate individual 20, 21, 22, 23, 24, 20-24 and 25-30 BED files from above raw and norm
##raw
awk -v OFS='\t' '{if ($8 == 20) print}' $BED_ROOT/${NAME}/${NAME}.bed > $BED_ROOT/${NAME}/${NAME}.20.bed
awk -v OFS='\t' '{if ($8 == 21) print}' $BED_ROOT/${NAME}/${NAME}.bed > $BED_ROOT/${NAME}/${NAME}.21.bed
awk -v OFS='\t' '{if ($8 == 22) print}' $BED_ROOT/${NAME}/${NAME}.bed > $BED_ROOT/${NAME}/${NAME}.22.bed
awk -v OFS='\t' '{if ($8 == 23) print}' $BED_ROOT/${NAME}/${NAME}.bed > $BED_ROOT/${NAME}/${NAME}.23.bed
awk -v OFS='\t' '{if ($8 == 24) print}' $BED_ROOT/${NAME}/${NAME}.bed > $BED_ROOT/${NAME}/${NAME}.24.bed
awk -v OFS='\t' '{if ($8 == 20 || $8 == 21 || $8 == 22) print}' $BED_ROOT/${NAME}/${NAME}.bed > $BED_ROOT/${NAME}/${NAME}.20to22.bed
awk -v OFS='\t' '{if ($8 == 23 || $8 == 24) print}' $BED_ROOT/${NAME}/${NAME}.bed > $BED_ROOT/${NAME}/${NAME}.23to24.bed
awk -v OFS='\t' '{if ($8 == 20 || $8 == 21 || $8 == 22 || $8 == 23 || $8 == 24) print}' $BED_ROOT/${NAME}/${NAME}.bed > $BED_ROOT/${NAME}/${NAME}.20to24.bed
awk -v OFS='\t' '{if ($8 == 25 || $8 == 26 || $8 == 27 || $8 == 28 || $8 == 29 || $8 == 30) print}' $BED_ROOT/${NAME}/${NAME}.bed > $BED_ROOT/${NAME}/${NAME}.25to30.bed
##norm
awk -v OFS='\t' '{if ($8 == 20) print}' $BED_ROOT/${NAME}/${NAME}.norm.bed > $BED_ROOT/${NAME}/${NAME}.norm.20.bed
awk -v OFS='\t' '{if ($8 == 21) print}' $BED_ROOT/${NAME}/${NAME}.norm.bed > $BED_ROOT/${NAME}/${NAME}.norm.21.bed
awk -v OFS='\t' '{if ($8 == 22) print}' $BED_ROOT/${NAME}/${NAME}.norm.bed > $BED_ROOT/${NAME}/${NAME}.norm.22.bed
awk -v OFS='\t' '{if ($8 == 23) print}' $BED_ROOT/${NAME}/${NAME}.norm.bed > $BED_ROOT/${NAME}/${NAME}.norm.23.bed
awk -v OFS='\t' '{if ($8 == 24) print}' $BED_ROOT/${NAME}/${NAME}.norm.bed > $BED_ROOT/${NAME}/${NAME}.norm.24.bed
awk -v OFS='\t' '{if ($8 == 20 || $8 == 21 || $8 == 22) print}' $BED_ROOT/${NAME}/${NAME}.norm.bed > $BED_ROOT/${NAME}/${NAME}.norm.20to22.bed
awk -v OFS='\t' '{if ($8 == 23 || $8 == 24) print}' $BED_ROOT/${NAME}/${NAME}.norm.bed > $BED_ROOT/${NAME}/${NAME}.norm.23to24.bed
awk -v OFS='\t' '{if ($8 == 20 || $8 == 21 || $8 == 22 || $8 == 23 || $8 == 24) print}' $BED_ROOT/${NAME}/${NAME}.norm.bed > $BED_ROOT/${NAME}/${NAME}.norm.20to24.bed
awk -v OFS='\t' '{if ($8 == 25 || $8 == 26 || $8 == 27 || $8 == 28 || $8 == 29 || $8 == 30) print}' $BED_ROOT/${NAME}/${NAME}.norm.bed > $BED_ROOT/${NAME}/${NAME}.norm.25to30.bed

#7. Make bedgraphs for normalized files for visualization (e.g. on IGV) (NOTE: will have to create a ${GENOME}.genome file listing the contigs and length of contigs in bp as two tab-delimited columns
##Compute scaling factor for normalization of bedgraphs; millions of genome-matching reads
SCALE=$(bc -l <<< 1/$MGMR)
echo $SCALE
###plus strand
awk -v OFS='\t' '{if ($8 == 20 || $8 == 21 || $8 == 22 || $8 == 23 || $8 == 24)print}' $BED_ROOT/${NAME}/${NAME}.sorted.bed | bedtools genomecov -i stdin -g $GENOME_DIR/${GENOME}.genome -bg -strand + -scale $SCALE > $BED_ROOT/${NAME}/${NAME}.20to24.norm.plus.bedgraph
awk -v OFS='\t' '{if ($8 == 20 || $8 == 21 || $8 == 22)print}' $BED_ROOT/${NAME}/${NAME}.sorted.bed | bedtools genomecov -i stdin -g $GENOME_DIR/${GENOME}.genome -bg -strand + -scale $SCALE > $BED_ROOT/${NAME}/${NAME}.20to22.norm.plus.bedgraph
awk -v OFS='\t' '{if ($8 == 23 || $8 == 24)print}' $BED_ROOT/${NAME}/${NAME}.sorted.bed | bedtools genomecov -i stdin -g $GENOME_DIR/${GENOME}.genome -bg -strand + -scale $SCALE > $BED_ROOT/${NAME}/${NAME}.23to24.norm.plus.bedgraph
awk -v OFS='\t' '{if ($8 == 25 || $8 == 26 || $8 == 27 || $8 == 28 || $8 == 29 || $8 == 30)print}' $BED_ROOT/${NAME}/${NAME}.sorted.bed | bedtools genomecov -i stdin -g $GENOME_DIR/${GENOME}.genome -bg -strand + -scale $SCALE > $BED_ROOT/${NAME}/${NAME}.25to30.norm.plus.bedgraph
##minus strand
awk -v OFS='\t' '{if ($8 == 20 || $8 == 21 || $8 == 22 || $8 == 23 || $8 == 24)print}' $BED_ROOT/${NAME}/${NAME}.sorted.bed | bedtools genomecov -i stdin -g $GENOME_DIR/${GENOME}.genome -bg -strand - -scale $SCALE > $BED_ROOT/${NAME}/${NAME}.20to24.norm.minus.bedgraph
awk -v OFS='\t' '{if ($8 == 20 || $8 == 21 || $8 ==22)print}' $BED_ROOT/${NAME}/${NAME}.sorted.bed | bedtools genomecov -i stdin -g $GENOME_DIR/${GENOME}.genome -bg -strand - -scale $SCALE > $BED_ROOT/${NAME}/${NAME}.20to22.norm.minus.bedgraph
awk -v OFS='\t' '{if ($8 == 23 || $8 == 24)print}' $BED_ROOT/${NAME}/${NAME}.sorted.bed | bedtools genomecov -i stdin -g $GENOME_DIR/${GENOME}.genome -bg -strand - -scale $SCALE > $BED_ROOT/${NAME}/${NAME}.23to24.norm.minus.bedgraph
awk -v OFS='\t' '{if ($8 == 25 || $8 == 26 || $8 == 27 || $8 == 28 || $8 == 29 || $8 == 30)print}' $BED_ROOT/${NAME}/${NAME}.sorted.bed | bedtools genomecov -i stdin -g $GENOME_DIR/${GENOME}.genome -bg -strand - -scale $SCALE > $BED_ROOT/${NAME}/${NAME}.25to30.norm.minus.bedgraph

#8. Quantify the number of reads overlapping annotated miRNA loci
module load BEDTools/2.26.0-foss-2016a

##Quantify 20, 21, 22, 23, 24 individually; 20-22, 23-24, 20-24 and 25-30 as groups
##raw
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.20.bed > $BED_ROOT/${NAME}/${NAME}.20.raw.on.$ANNO.bed
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.21.bed > $BED_ROOT/${NAME}/${NAME}.21.raw.on.$ANNO.bed
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.22.bed > $BED_ROOT/${NAME}/${NAME}.22.raw.on.$ANNO.bed
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.23.bed > $BED_ROOT/${NAME}/${NAME}.23.raw.on.$ANNO.bed
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.24.bed > $BED_ROOT/${NAME}/${NAME}.24.raw.on.$ANNO.bed
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.20to22.bed > $BED_ROOT/${NAME}/${NAME}.20to22.raw.on.$ANNO.bed
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.23to24.bed > $BED_ROOT/${NAME}/${NAME}.23to24.raw.on.$ANNO.bed
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.20to24.bed > $BED_ROOT/${NAME}/${NAME}.20to24.raw.on.$ANNO.bed
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.25to30.bed > $BED_ROOT/${NAME}/${NAME}.25to30.raw.on.$ANNO.bed
##norm
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.norm.20.bed > $BED_ROOT/${NAME}/${NAME}.20.norm.on.$ANNO.bed
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.norm.21.bed > $BED_ROOT/${NAME}/${NAME}.21.norm.on.$ANNO.bed
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.norm.22.bed > $BED_ROOT/${NAME}/${NAME}.22.norm.on.$ANNO.bed
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.norm.23.bed > $BED_ROOT/${NAME}/${NAME}.23.norm.on.$ANNO.bed
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.norm.24.bed > $BED_ROOT/${NAME}/${NAME}.24.norm.on.$ANNO.bed
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.norm.20to22.bed > $BED_ROOT/${NAME}/${NAME}.20to22.norm.on.$ANNO.bed
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.norm.23to24.bed > $BED_ROOT/${NAME}/${NAME}.23to24.norm.on.$ANNO.bed
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.norm.20to24.bed > $BED_ROOT/${NAME}/${NAME}.20to24.norm.on.$ANNO.bed
bedtools map -s -F $FRAG_OVERLAP -c 5 -o sum -null 0 -a $ANNO_ROOT/$ANNO.bed -b $BED_ROOT/${NAME}/${NAME}.norm.25to30.bed > $BED_ROOT/${NAME}/${NAME}.25to30.norm.on.$ANNO.bed

##combine quantifications BED files into one document: Name, 20, 21, 22, 23, 24, 20-24, 25-30

echo -e 'name\t20.raw\t21.raw\t22.raw\t23.raw\t24.raw\t20to22.raw\t23to24.raw\t20to24.raw\t25to30.raw\t20.norm\t21.norm\t22.norm\t23.norm\t24.norm\t20to22.norm\t23to24.norm\t20to24.norm\t25to30.norm' > $BED_ROOT/${NAME}/${NAME}.quant.on.$ANNO.tsv

#generate tmp files for later appending
awk -v OFS='\t' '{print $4}' $ANNO_ROOT/$ANNO.bed > $BED_ROOT/${NAME}/$ANNO.bed.tmp
#raw
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.20.raw.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.20.raw.on.$ANNO.bed.tmp
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.21.raw.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.21.raw.on.$ANNO.bed.tmp
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.22.raw.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.22.raw.on.$ANNO.bed.tmp
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.23.raw.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.23.raw.on.$ANNO.bed.tmp
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.24.raw.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.24.raw.on.$ANNO.bed.tmp
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.20to22.raw.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.20to22.raw.on.$ANNO.bed.tmp
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.23to24.raw.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.23to24.raw.on.$ANNO.bed.tmp
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.20to24.raw.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.20to24.raw.on.$ANNO.bed.tmp
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.25to30.raw.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.25to30.raw.on.$ANNO.bed.tmp
#norm
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.20.norm.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.20.norm.on.$ANNO.bed.tmp
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.21.norm.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.21.norm.on.$ANNO.bed.tmp
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.22.norm.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.22.norm.on.$ANNO.bed.tmp
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.23.norm.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.23.norm.on.$ANNO.bed.tmp
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.24.norm.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.24.norm.on.$ANNO.bed.tmp
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.20to22.norm.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.20to22.norm.on.$ANNO.bed.tmp
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.23to24.norm.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.23to24.norm.on.$ANNO.bed.tmp
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.20to24.norm.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.20to24.norm.on.$ANNO.bed.tmp
awk -v OFS='\t' '{print $7}' $BED_ROOT/${NAME}/${NAME}.25to30.norm.on.$ANNO.bed > $BED_ROOT/${NAME}/${NAME}.25to30.norm.on.$ANNO.bed.tmp

paste $BED_ROOT/${NAME}/$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.20.raw.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.21.raw.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.22.raw.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.23.raw.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.24.raw.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.20to22.raw.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.23to24.raw.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.20to24.raw.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.25to30.raw.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.20.norm.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.21.norm.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.22.norm.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.23.norm.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.24.norm.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.20to22.norm.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.23to24.norm.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.20to24.norm.on.$ANNO.bed.tmp \
$BED_ROOT/${NAME}/${NAME}.25to30.norm.on.$ANNO.bed.tmp \
>> $BED_ROOT/${NAME}/${NAME}.quant.on.$ANNO.tsv

#remove tmp files
rm $BED_ROOT/${NAME}/*tmp

echo 'sRNA pipe complete!'
