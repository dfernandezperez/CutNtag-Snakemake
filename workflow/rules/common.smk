def is_single_end(sample):
    return pd.isnull(units.loc[(sample), "fq2"][0])

def is_spike(sample):
    return SAMPLES.SPIKE[sample] == True

def has_input(sample):
    return SAMPLES.INPUT[sample] != "NoInput"


# Get raw or trimmed reads based on trimming configuration
def get_fq(wildcards):
    if config["trimming"]:
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("{tmp}/fastq/trimmed/{sample}.{group}.fastq.gz", group=[1, 2], **wildcards, tmp = config["tmp"])
        # single end sample
        return "{tmp}/fastq/trimmed/{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)
    else:
        # no trimming, use raw reads
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("{tmp}/fastq/{sample}.{group}.fastq.gz", group=[1, 2], **wildcards, tmp = config["tmp"])
        # single end sample
        return "{tmp}/fastq/{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)


def get_fq_spike(wildcards):
    if is_spike(**wildcards):
        if config["trimming"]:
            if not is_single_end(**wildcards):
                # paired-end sample
                return expand("{tmp}/fastq/trimmed{sample}.{group}.fastq.gz", group=[1, 2], **wildcards, tmp = config["tmp"])
            # single end sample
            return "{tmp}/fastq/trimmed{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)
        else:
            # no trimming, use raw reads
            if not is_single_end(**wildcards):
                # paired-end sample
                return expand("{tmp}/fastq/{sample}.{group}.fastq.gz", group=[1, 2], **wildcards, tmp = config["tmp"])
            # single end sample
            return "{tmp}/fastq/{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)


# Get raw or trimmed reads based on trimming configuration. Used for fastqc
def get_fq_forward(wildcards):
    if config["trimming"]:
        if not is_single_end(**wildcards):
            # paired-end sample
            return "{tmp}/fastq/trimmed/{sample}.1.fastq.gz".format(**wildcards, tmp = config["tmp"])
        # single end sample
        return "{tmp}/fastq/trimmed/{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)
    else:
        # no trimming, use raw reads
        if not is_single_end(**wildcards):
            # paired-end sample
            return "{tmp}/fastq/{sample}.1.fastq.gz".format(**wildcards, tmp = config["tmp"])
        # single end sample
        return "{tmp}/fastq/{sample}.se.fastq.gz".format(tmp = config["tmp"], **wildcards)


def get_bam(wildcards):
    if not is_spike(**wildcards):
        return "results/02aln/{sample}.bam.tmp".format(**wildcards)
    return "results/02aln/{sample}.bam.tmp.clean".format(**wildcards)


def set_read_extension(wildcards):
    if is_single_end(wildcards.sample):
        return "--extendReads " + str(config['bam2bigwig']['read_extension'])
    return "--extendReads"


def get_bam_spike(wildcards):
    if not is_spike(wildcards.sample):
        return { "case": "results/02aln/{sample}.bam".format(sample=wildcards.sample) }
    else :
        return { "case": "results/02aln/{sample}.bam".format(sample=wildcards.sample),
             "spike": "results/02aln_dm/{sample}_spike.bam.clean".format(sample=wildcards.sample) }


def set_reads(wildcards, input):
        n = len(input)
        if n == 1:
            reads = "{}".format(*input)
            return reads
        else:
            reads = config["params"]["bowtie2"]["pe"] + " -1 {} -2 {}".format(*input)
            return reads


def set_reads_spike(wildcards, input):
        n = len(input)
        assert n == 1 or n == 2, "input->sample must have 1 (sample) or 2 (sample + spike) elements"
        if n == 1:
            reads = "workflow/scripts/bam2bigwig.py"
            return reads
        if n == 2:
            reads = "workflow/scripts/bam2bigwig_spike.py --spike {}".format(input.spike)
            return reads