import glob
from os.path import join

ref_dir = "/home/acheng/data/reference/"
STAR_INDEX = ref_dir + "hg19/star/refseq"

S_STAR = "/home/acheng/software/STAR-2.5.2a/bin/Linux_x86_64/STAR"
S_SAMTOOLS = "/home/acheng/software/samtools-1.1/samtools"
S_PYTHON = "/home/acheng/software/miniconda3/envs/karrseq/bin/python"
S_SEQPREP = "/home/acheng/software/SeqPrep/SeqPrep"
S_SORT = "/usr/bin/sort"
S_PAIRIX = "/home/acheng/software/pairix/bin/pairix"
S_BGZIP = "/home/acheng/software/pairix/bin/bgzip"

# =================
# Targets
# =================
SAMPLES = [ sample.split("/")[-1].replace("_P1.fastq.gz", "")
	    for sample in glob.glob("data/fastq/*_P1.fastq.gz") ]

rule targets_all:
    input:
        expand("data/inner/{genome}/MAPQ{mapq}_SPAN{span}/{sample}.{todedup}.pairs.gz",
               sample=SAMPLES,
               mapq=["1"],
               span=["0"],
               todedup=["dedup"],
               genome=["hrefseq"])


# =====================
# Rules for Single-end (STAR)
# =====================
rule merge_fastq:
    input:
        read1 = "data/fastq/{sample}_P1.fastq.gz",
        read2 = "data/fastq/{sample}_P2.fastq.gz"
    output:
        read1 = "data/fastq/merge/{sample}_unmerge_P1.fastq.gz",
        read2 = "data/fastq/merge/{sample}_unmerge_P2.fastq.gz",
        read3 = "data/fastq/merge/{sample}_reject_P1.fastq.gz",
        read4 = "data/fastq/merge/{sample}_reject_P2.fastq.gz",
        merge = "data/fastq/merge/{sample}_merge.fastq.gz"
    shell:
        """
	{S_SEQPREP} \
            -f {input.read1} \
            -r {input.read2} \
            -1 {output.read1} \
            -2 {output.read2} \
            -3 {output.read3} \
            -4 {output.read4} \
            -s {output.merge}
        """


rule align_with_STAR:
    input: "data/fastq/merge/{sample}_merge.fastq.gz"
    output:
        log = "data/bam/{genome}/{sample}.log",
	aligned = "data/bam/{genome}/{sample}_Aligned.sortedByCoord.out.bam",
        chimeric = "data/bam/{genome}/{sample}_Chimeric.out.sam"
    params:
        prefix = "data/bam/{genome}/{sample}",
        index = STAR_INDEX
    threads: 4
    shell:
        """
        {S_STAR} \
            --runMode alignReads \
            --runThreadN {threads} \
            --genomeDir {params.index} \
            --readFilesIn {input} \
            --readFilesCommand zcat \
            --outFileNamePrefix {params.prefix}_ \
            --outReadsUnmapped Fastq \
            --outFilterMultimapNmax 100 \
            --outSAMtype BAM SortedByCoordinate \
            --outSAMattributes All \
            --alignIntronMin 1 \
            --scoreGapNoncan -4 \
            --scoreGapATAC -4 \
            --chimSegmentMin 15 \
            --chimJunctionOverhangMin 15 \
            --limitOutSJcollapsed 10000000 \
            --limitIObufferSize 1500000000
        {S_SAMTOOLS} index {output.aligned}
        touch {output.log}
        """

rule star_single_end_chimeric_reads_to_pairs:
    input:
        log = "data/bam/{genome}/{sample}.log",
        aligned = "data/bam/{genome}/{sample}_Aligned.sortedByCoord.out.bam",
        chimeric = "data/bam/{genome}/{sample}_Chimeric.out.sam"
    output: "data/pairs/{genome}/MAPQ{mapq}_SPAN{span}/{sample}.txt.gz"
    params:
        mapq = "{mapq}",
        span = "{span}",
        tmp1 = "data/pairs/{genome}/MAPQ{mapq}_SPAN{span}/{sample}.tmp1",
        tmp2 = "data/pairs/{genome}/MAPQ{mapq}_SPAN{span}/{sample}.tmp2"
    threads: 4
    shell:
        """
        {S_SAMTOOLS} view {input.aligned} \
            | {S_PYTHON} src/get_STAR_reads.py Aligned {params.mapq} {params.span} > {params.tmp1}
        cat {input.chimeric} \
            | {S_PYTHON} src/get_STAR_reads.py Chimeric {params.mapq} {params.span} > {params.tmp2}
        cat {params.tmp1} {params.tmp2} \
            | {S_SORT} -k2,2 -k4,4 -k3,3n -k5,5n -k10,10n -k11,11n --parallel={threads} \
            | gzip -c > {output}
        rm {params.tmp1} {params.tmp2}
        """

rule star_single_end_remove_duplicates:
    input: "data/pairs/{genome}/MAPQ{mapq}_SPAN{span}/{sample}.txt.gz"
    output:
        dedup = "data/pairs/{genome}/MAPQ{mapq}_SPAN{span}/{sample}.dedup.txt.gz",
        bed = "data/pairs/{genome}/MAPQ{mapq}_SPAN{span}/{sample}.dedup.bed",
        pairs = "data/pairs/{genome}/MAPQ{mapq}_SPAN{span}/{sample}.dedup.pairs.gz"
    params:
        mapq = "{mapq}",
        todedup = "dedup"
    shell:
        """
        zcat {input} | {S_PYTHON} src/remove_duplicates.py {params.todedup} | gzip -c > {output.dedup}
        zcat {output.dedup} | {S_PYTHON} src/pairs_to_bed.py {params.mapq} > {output.bed}
        zcat {output.dedup} | cut -f1-7 | {S_BGZIP} -c > {output.pairs}
        sleep 120
        {S_PAIRIX} {output.pairs}
        """

rule star_single_end_ligation:
    input: "data/pairs/{genome}/MAPQ{mapq}_SPAN{span}/{sample}.dedup.txt.gz"
    output: "data/inner/{genome}/MAPQ{mapq}_SPAN{span}/{sample}.dedup.pairs.gz"
    threads: 4
    shell:
        """
        zcat {input} \
            | awk '{{print $1"\t"$2"\t"$3+$10"\t"$4"\t"$5"\t"$6"\t"$7}}' \
            | awk '{{if($2==$4){{if($3>$5){{print $1"\t"$4"\t"$5"\t"$2"\t"$3"\t"$7"\t"$6}}else{{print $0}}}}else{{print $0}}}}' \
            | {S_SORT} -k2,2 -k4,4 -k3,3n -k5,5n --parallel={threads} \
            | {S_BGZIP} -c > {output}
        sleep 60
        {S_PAIRIX} {output}
        """
