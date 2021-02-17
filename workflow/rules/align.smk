rule align:
    input:
        get_fq
    output:
         bam   = temp("results/02aln/{sample}.bam.tmp"),
         index = temp("results/02aln/{sample}.bam.tmp.bai")
    threads:
        CLUSTER["align"]["cpu"]
    params:
        index  	     = config["ref"]["index"],
        bowtie2 	 = config["params"]["bowtie2"]["global"],
        samblaster   = config["params"]["samblaster"],
        reads  	     = set_reads,
        samtools_mem = config["params"]["samtools"]["memory"]
    message:
        "Aligning {input} with parameters {params.bowtie2}"
    log:
       align   = "results/00log/alignments/{sample}.log",
       rm_dups = "results/00log/alignments/rm_dup/{sample}.log",
    benchmark:
        "results/.benchmarks/{sample}.align.benchmark.txt"
    shell:
        """
        bowtie2 -p {threads} {params.bowtie2} -x {params.index} {params.reads} 2> {log.align} \
        | samblaster {params.samblaster} 2> {log.rm_dups} \
        | samtools view -q2 -Sb -F 4 - \
        | samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
        samtools index {output.bam}
        """


rule align_spike:
    input:
        get_fq_spike
    output:
        bam   = temp("results/02aln_dm/{sample}_spike.bam"),
        index = temp("results/02aln_dm/{sample}_spike.bam.bai")
    threads:
        CLUSTER["align"]["cpu"]
    params:
        index        = config["ref"]["index_spike"],
        bowtie2      = config["params"]["bowtie2"]["global"],
        samblaster   = config["params"]["samblaster"],
        reads        = set_reads,
        samtools_mem = config["params"]["samtools"]["memory"]
    message:
        "Aligning {input} with parameters {params.bowtie2}"
    log:
       align   = "results/00log/alignments/{sample}_spike.log",
       rm_dups = "results/00log/alignments/rm_dup/{sample}_spike.log",
    benchmark:
        "results/.benchmarks/{sample}.alignSpike.benchmark.txt"
    shell:
        """
        bowtie2 -p {threads} {params.bowtie2} -x {params.index} {params.reads} 2> {log.align} \
        | samblaster {params.samblaster} 2> {log.rm_dups} \
        | samtools view -Sb -F 4 - \
        | samtools sort -m {params.samtools_mem}G -@ {threads} -T {output.bam}.tmp -o {output.bam} - 2>> {log.align}
        samtools index {output.bam}
        """


rule clean_spike:
    input:
        mm          = "results/02aln/{sample}.bam.tmp",
        spike       = "results/02aln_dm/{sample}_spike.bam",
        mm_index    = "results/02aln/{sample}.bam.tmp.bai",
        spike_index = "results/02aln_dm/{sample}_spike.bam.bai",
    output:
        mm    = temp("results/02aln/{sample}.bam.tmp.clean"),
        spike = "results/02aln_dm/{sample}_spike.bam.clean"
    log:
        "results/00log/alignments/{sample}.removeSpikeDups"
    shell:
        """
        python workflow/scripts/remove_spikeDups.py {input} &> {log}      
        mv {input.mm}.temporary {output.mm}; mv {input.spike}.temporary {output.spike}
        samtools index {output.spike}
        """

# Dummy rule to change the name of the bam files to be able to 
# have the same name structure in spike-in and non-spiked samples
rule update_bam:
    input:
        get_bam
    output:
        "results/02aln/{sample}.bam",
    log:
        "results/00log/alignments/{sample}.update_bam"
    shell:
        """
        cp {input} {output}
        samtools index {output} 2>> {log}
        """

rule bam2bigwig:
    input: 
        unpack(get_bam_spike)
    output:  
        "results/06bigwig/{sample}.bw"
    params: 
        read_exten = set_read_extension,
        reads      = set_reads_spike,
        params     = config["bam2bigwig"]["other"]
    log: 
        "results/00log/bam2bw/{sample}_bigwig.bam2bw"
    threads: 
        CLUSTER["bam2bigwig"]["cpu"]
    message: 
        "making bigwig for sample {wildcards.sample}"
    shell:
        """
        python {params.reads} \
        --case {input.case} \
        --output {output} \
        --threads {threads} \
        --otherParams {params.read_exten} {params.params} &> {log}
        """


rule bigwig2server:
    input: 
        bw         = "results/06bigwig/{sample}.bw",
        samblaster = "results/00log/alignments/rm_dup/{sample}.log",
        bowtie     = "results/00log/alignments/{sample}.log"
    output:
        temp("results/temp_file_{sample}_{control}.txt")
    params:
        user     = lambda wildcards : SAMPLES.USER[wildcards.sample],
        antibody = lambda wildcards : SAMPLES.AB[wildcards.sample],
        genome   = lambda wildcards : SAMPLES.GENOME[wildcards.sample],
        run      = lambda wildcards : SAMPLES.RUN[wildcards.sample],
        chip     = lambda wildcards : str("CutNtag") if SAMPLES.SPIKE[wildcards.sample] == False else str("CutNtagSpike")
    run:
        # Get number of removed reported reads by bowtie
        with open(input.bowtie,"r") as fi:
            for ln in fi:
                lines = str.split(ln)
                if "exactly" in lines:
                    nreads = int(lines[0])
                if ">1" in lines:
                    nreads += int(lines[0])
        # Get number of removed reads
        with open(input.samblaster,"r") as fi:
            for ln in fi:
                if ln.startswith("samblaster: Marked "):
                    removed_reads = int( str.split(ln)[2] )

        # Total number of final reads is reported by bowtie minus duplicated removed
        total_reads = nreads-removed_reads

        shell(
            "cp {input} \
            /hpcnfs/data/DP/UCSC_tracks/Data/bigWig/{sample}_{control}_{user}_{nreads}_{chip}_{antibody}_{genome}_{run}.bigWig".format(
            input    = input.bw,
            sample   = wildcards.sample,
            control  = wildcards.control,
            user     = params.user,
            nreads   = total_reads,
            chip     = params.chip,
            antibody = params.antibody,
            genome   = params.genome,
            run      = params.run)
            )
        shell("touch {output}".format(output = output))
