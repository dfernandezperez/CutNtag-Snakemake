rule bam2bedgraph:
    input: 
        unpack(get_bam_spike)
    output:  
        temp("results/02aln/bedgraph/{sample}.bedgraph")
    params: 
        read_exten = set_read_extension,
        reads      = set_reads_spike,
        params     = config["bam2bigwig"]["other"] # Same params as bam2bw
    log: 
        "results/00log/bam2bedgraph/{sample}.log"
    threads: 
        CLUSTER["bam2bedgraph"]["cpu"]
    message: 
        "making bedgraph for sample {wildcards.sample}"
    shell:
        """
        python {params.reads} \
        --case {input.case} \
        --output {output} \
        --threads {threads} \
        --otherParams --outFileFormat bedgraph {params.read_exten} {params.params} &> {log}
        """

rule seacr_noIgG:
    input:
        "results/02aln/bedgraph/{sample}.bedgraph"
    output:
        peaks = "results/03seacr/{sample}/{sample}_{filt}.peaks.stringent.bed"
    params:
        mode       = config["params"]["seacr"]["mode"],
        out_prefix = "results/03seacr/{sample}/{sample}_{filt}.peaks"
    log:
        "results/00log/seacr/{sample}_top{filt}.log"
    shadow:
        "minimal"
    shell:
        """
        bash /SEACR_1.3.sh {input} {wildcards.filt} non {params.mode} {params.out_prefix} 2> {log}
        """