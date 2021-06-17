import sys

FILES, = glob_wildcards("{file}.fastq")
if len(FILES) == 0:
    sys.exit()

BARCODES = config["barcodes"]

REFERENCE = "/mnt/StudentFiles/2020-21/Project5/data/influenza/influenza_a_hongkong_h9n2.fasta"
REFERENCE_GFF = "/mnt/StudentFiles/2020-21/Project5/data/covid/NC_045512.2_edited.gff3"
HG38 = "/mnt/StudentFiles/2020-21/Project5/data/GRCh38_latest_genomic.fna"
COVERAGE = ["3"]

rule all:
    input:
        expand("blast/blast_on_barcode{bc}_consensus_{coverage}x.txt", bc=BARCODES, coverage=COVERAGE),
	expand("nanoplot/barcode{bc}", bc=BARCODES)
        
rule demux_files:
    priority: 8
    input:
        "{file}.fastq"
    output:
        barcode = expand("separate/{{file}}/BC{bc}.fastq", bc=BARCODES)
    params:
        barcodes = expand("BC{bc}_rapid", bc=BARCODES)
    threads: 1
    shell:
        """
        /mnt/StudentFiles/2020-21/Project5/Porechop1/porechop-runner.py -i {input} -b separate/{wildcards.file} \
        --format fastq --require_two_barcodes --discard_unassigned -t {threads} -v 0 -l {params.barcodes}
        touch {output.barcode}
        """

rule NanoPlot:
    priority: 7
    input:
        "demultiplexed/barcode{bc}.fastq"
    output:
        "nanoplot/barcode{bc}"
    threads: 1
    shell:
        """
        NanoPlot --fastq {input} -o {output} -t {threads}
        """

rule merge_barcodes:
    priority: 6
    input:
        expand("separate/{file}/BC{{bc}}.fastq", file=FILES)
    output:
        "demultiplexed/barcode{bc,\d+}.fastq"
    threads: 1
    shell:
        """
        cat {input} > {output}
        """
        
rule cut_adapters:
    priority: 5
    input:
        "demultiplexed/barcode{bc}.fastq"
    output:
        "trimmed/barcode{bc}_trimmed.fastq"
    threads: 1
    shell:
        """
        cutadapt -u 30 -u -30 -o {output} {input} -m 75 -j {threads} --quiet
        """

rule filter_human_out:
    priority: 4
    input:
        trimmed_fastq="trimmed/barcode{bc}_trimmed.fastq",
        reference=HG38
    output:
        "filtered/barcode{bc,\d+}_filtered.fastq"
    threads: 8
    shell:
        """
        minimap2 -Y -t {threads} -x map-ont -a {input.reference} {input.trimmed_fastq} 2> /dev/null | samtools fastq -f 4 - 2> /dev/null > {output}
        """

rule map_to_reference:
    priority: 3
    input:
        filtered_fastq="filtered/barcode{bc}_filtered.fastq",
        reference=REFERENCE
    output:
        "mapped/barcode{bc,\d+}_mapped.bam"
    threads: 4
    shell:
        """
        minimap2 -Y -t {threads} -x map-ont -a {input.reference} {input.filtered_fastq} 2> /dev/null | samtools view -bF 4 - | samtools sort -@ {threads} - > {output}
        """

rule create_consensus:
    priority: 2
    input:
        "mapped/barcode{bc}_mapped.bam"
    output:
        "consensus/barcode{bc,\d+}_consensus_{coverage,\d+}x.fasta"
    threads: 1
    shell:
        """
        samtools index -@ {threads} {input}
        /mnt/StudentFiles/2020-21/Project5/covid_snake/helper_scripts/bam2consensus.py -i {input} -o {output} -d {wildcards.coverage} -g 1
        """

rule blast:
    priority: 1
    input:
        "consensus/barcode{bc,\d+}_consensus_{coverage}x.fasta"
    output:
        "blast/blast_on_barcode{bc,\d+}_consensus_{coverage}x.csv"
    threads: 1
    shell:
        """
        blastn -query {input} -db /mnt/StudentFiles/2020-21/Project5/data/influenza/database/influenza -out {output} -outfmt "10 qseqid stitle evalue qcovs" -evalue 0.00001 -max_target_seqs 1
        """


