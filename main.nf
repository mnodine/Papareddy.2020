#!/usr/bin/env nextflow


/**********************
* parameters
**********************/
params.output = "results"
params.files = "/PATH/TO/FASTQs/*.fastq"
params.fasta_dna = "./Plotnikova.2019/Arabidopsis_thaliana.TAIR10.dna.toplevel.fa"
params.adapt = "AGATCGGAAGA"
//STAR
params.mismatch = 0
params.filterMNmax = 100
params.intronmax = 1
params.aligntype = "EndToEnd" 
params.outputtype = "BAM SortedByCoordinate"
params.unmapped = "Fastx"
params.flag = "AllBestScore"
params.trim=true
/**********************
* set up channels
********************/
files = Channel.fromPath(params.files).map { file -> [ id:file.baseName,file:file] }
.ifEmpty { error "Cannot find any reads matching: ${params.files}" }

fasta_dna = file(params.fasta_dna)
if( !fasta_dna.exists() ) exit 1, "Missing dna fasta file: ${fasta_dna}"



files.into{to_trim;for_fastQC}



/**********************
* fastqc
********************/
process fastqc {
    tag "$name"
    publishDir "${params.output}/fastqc", mode: 'copy',
        saveAs: {filename -> filename.indexOf(".zip") > 0 ? "zips/$filename" : "$filename"}

    input:
    set val(name), file(fq) from for_fastQC

    output:
    file '*_fastqc.{zip,html}' into fastqc_results

    script:
    """
    fastqc -q $fq
    """
}


/********************
* trim adapter
* use only 18-30 bp long fragments
* see https://github.com/Gregor-Mendel-Institute/Plotnikova.2019/blob/master/sRNA.align.git.sh
*********************/

if (!params.trim ){
	to_trim.set {trimmed}
} else {

process trim {
publishDir "$params.output/trimmed", mode: 'copy'
input:
set id, file(fq) from to_trim

output:
set id, file("${id}_trimmed.fq") into trimmed
 file("${id}.txt") into trimmed_log
script:
"""
cutadapt -a ${params.adapt}  -m 18 -M 30  --trimmed-only ${fq} >  ${id}_trimmed.fq 2> ${id}.txt

"""
}
}
trimmed.into{ tobowtie; tostar}

/*********************
* STAR INDEX
*********************/

process star_index {

  input:
  file fasta from fasta_dna

  output:
  file "star_index" into index

  script:
  """
  mkdir -p star_index
  STAR --runThreadN 8 --runMode genomeGenerate --genomeDir star_index --genomeFastaFiles ${fasta}
  """
}


/********************
* STAR
* similar to https://github.com/Gregor-Mendel-Institute/Plotnikova.2019/blob/master/sRNA.align.git.sh
*********************/

process star_align {
publishDir "$params.output/alignments", mode: 'copy'
  input:
  file index from index
  set id, file(fq) from tostar

  output:
  set id, file("${id}Aligned.sortedByCoord.out.bam") into star_bam
  set id, file("${id}Log.final.out") into star_log

  script:
  """
  STAR --genomeDir $index \
  --readFilesIn $fq \
  --outFileNamePrefix ${id} \
  --outFilterMismatchNmax $params.mismatch \
  --outFilterMultimapNmax $params.filterMNmax \
  --outSAMtype $params.outputtype \
  --outReadsUnmapped $params.unmapped \
  --alignIntronMax $params.intronmax \
  --alignEndsType $params.aligntype \
  --outSAMprimaryFlag $params.flag
  """
}


/**************************
* bamtosam
***************************/
process bamtosam {

input:
  
  set id, file(bam) from star_bam

  output:
  set id, file("${id}.sam") into bam2sam
  
script:
  """
  samtools view -h -o ${id}.sam ${bam}
  
  """
 
}

bam2sam.into{split_sam;for_readmapIO}


size = [20,21,22,23,24]

/**************************
* split by sam and bam files by length

***************************/

process sam_by_size {
publishDir "$params.output/sbam_by_length", mode: 'copy'

input:


  set id, file(sam) from split_sam
  each ntlength from size

 output:
  set id, file("${id}.${ntlength}.bam") into sam_sized

shell:
'''
awk '/^@/ || length($10) == !{ntlength}' !{sam}  | samtools view -S -b -h -o !{id}.!{ntlength}.bam
'''

}





/*****************************
* BAM 2 BW
****************************/

process bam2bw {
publishDir "$params.output/bam_bw", mode: 'copy'


input:
set var, file(bam) from sam_sized


output:
set var, file("${bam}.bw")

script:
"""
export TMPDIR=\$(pwd)
samtools index -b ${bam}
bamCoverage -b ${bam} -o ${bam}.bw  --normalizeUsing 'RPGC'  --effectiveGenomeSize 119146348  -bs 25 --smoothLength 50
"""
}




/**************************
* readmapIO
***************************/
process readmapIO {
publishDir "$params.output/mapIO", mode: 'copy'

input:
  
  set id, file(sam) from for_readmapIO

output:
  set id, file("${id}.BODY.temp.bed") into rescuedIO_bed
  set id, file("${id}_BODY_minus.bedgraph") into rescuedIO_BODY_minus_bedgraph
  set id, file("${id}_BODY_minus_coverage.bedgraph") into rescuedIO_BODY_minus_coverage_bedgraph
  set id, file("${id}_BODY_minus_untemp.bedgraph") into rescuedIO_minus_untemp_bedgraph
  set id, file("${id}_BODY_plus.bedgraph") into rescuedIO_BODY_plus_bedgraph
  set id, file("${id}_BODY_plus_coverage.bedgraph") into rescuedIO_BODY_plus_coverage_bedgraph
  set id, file("${id}_BODY_plus_untemp.bedgraph") into rescuedIO_BODY_plus_untemp_bedgraph
  

script:
  """
  python $baseDir/Plotnikova.2019/readmapIO.py -S ${id} -I ${id}.sam --bed -F $fasta_dna  --allow_naive --untemp_out ACGTN 

  
  """
  
}

/**************************
* bedtoolsort
***************************/
process bedSort {


input:
    set id, file(bed) from rescuedIO_bed

output:
  set id, file("${id}.sorted.bed") into sortedBed
  
script:
  """
 bedtools sort -i ${id}.BODY.temp.bed > ${id}.sorted.bed
  
  """
  
}



/**************************
* collapse overlapping reads
***************************/
process bedCollaps {

publishDir "$params.output/collapesedreads", mode: 'copy'
input:
  
  set id, file(bed) from sortedBed

  output:
  set id, file("${id}.collapsed.bed") into collapsedBed
  
script:
  """
  python $baseDir/Plotnikova.2019/bed_collapse.py ${id}.sorted.bed > ${id}.collapsed.bed
  
  """
  
}




/**************************
* hit normalization
*get total number of genome-mapping reads from BG file and use it to normalize BG hit-normalized reads
***************************/

process hitnormalize {

publishDir "$params.output/hitnorm", mode: 'copy'
input:
  
  set id, file(bed) from collapsedBed

 output:
  set id, file("${id}*.bed") into hitnormalizedBed

shell:
'''
    GMR=$(awk -v OFS='\t' '{sum+=$5}END{print sum}' !{bed})
    MGMR=$(awk '{print +$1 }' <<< $GMR)
    MGMR=$(bc -l <<< $MGMR/1000000)
   
    awk -v OFS='\t' -v MGMR="$MGMR" '{print $1,$2,$3,$4,$5/MGMR,$6,$7,$8}' !{bed} > !{id}.norm.bed
'''

}


size = [20,21,22,23,24]

/**************************
* split by read length
*small RNA size class noirmalide RPM
***************************/

process sizeclass {

publishDir "$params.output/sizeclass", mode: 'copy'
input:
  
  set id, file(bed) from hitnormalizedBed
  each length from size
 output:
  set id, file("${id}.norm.${length}.bed") into size_class

shell:
'''
awk -v OFS='\t' '{if ($8 == !{length}) print}' !{bed} > !{id}.norm.!{length}.bed
'''

}


/**************************
*TE derived RPM of siRNAs by size
***************************/

process TEmap {
publishDir "$params.output/TEmapped", mode: 'copy'

input:

  set id, file(bed) from size_class

  output:
  set id, file("${id}.TEmapped") into TEmapped

script:
  """
  bedtools map -a $baseDir/Plotnikova.2019/TEs.siRNA.clusterID.sorted -b ${bed} -o sum -g $baseDir/Plotnikova.2019/ara.genome -null 0  > ${id}.TEmapped 

  """

}














