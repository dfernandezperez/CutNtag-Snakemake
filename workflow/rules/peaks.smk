def input_macs2(wildcards):
    if has_input(wildcards.sample):
        return { "case": "results/02aln/{sample}.bam".format(sample=wildcards.sample),
             "reference": "results/02aln/{control}.bam".format(control=wildcards.control) }
    # Set reference as the case because otherwise I can't use input.reference inside shell, it complains that reference doesn't exist if the sample has no input.
    # Since in case there's not input that macs2 code is not going to be executed, is not a problem to set this dummy variable as a workaround.
    else:
        return { "case": "results/02aln/{sample}.bam".format(sample=wildcards.sample),
            "reference": "results/02aln/{sample}.bam".format(sample=wildcards.sample) }


rule call_peaks:
    input: 
        unpack(input_macs2)
    output: 
        narrowPeak = "results/03peak_macs2/{sample}_{control}/{sample}_peaks.narrowPeak",
        xls        = "results/03peak_macs2/{sample}_{control}/{sample}_peaks.xls"
    log:
        "results/00log/macs2/{sample}_{control}_macs2.log"
    params:
        out_dir      = "results/03peak_macs2/{sample}_{control}",
        macs2_params = config["params"]["macs2"]["pk_calling"],
        pvalue       = config["params"]["macs2"]["pvalue"],
        gsize        = config["params"]["macs2"]["gsize"],
        paired_end   = lambda w: "--format BAM --nomodel" if is_single_end(w.sample) else "--format BAMPE"
    message: 
        "call_peaks macs2 for sample {input.case}"
    shell:
        """
        if [[ {wildcards.control} != "NoInput" ]]
        then
            macs2 callpeak {params.paired_end} \
                --treatment {input.case} \
                --control {input.reference} \
                --gsize {params.gsize} \
                --outdir {params.out_dir} \
                --name {wildcards.sample} \
                --pvalue {params.pvalue} \
                {params.macs2_params} 2> {log}
        else
            macs2 callpeak {params.paired_end} \
                --treatment {input.case} \
                --gsize {params.gsize} \
                --outdir {params.out_dir} \
                --name {wildcards.sample} \
                --pvalue {params.pvalue} \
                {params.macs2_params} 2> {log}
        fi                   
        """

rule call_peaks_broad:
    input: 
        unpack(input_macs2)
    output: 
        broadPeak  = "results/03peak_macs2/{sample}_{control}/broad/{sample}_peaks.broadPeak",
        gappedPeak = "results/03peak_macs2/{sample}_{control}/broad/{sample}_peaks.gappedPeak",
        xls        = "results/03peak_macs2/{sample}_{control}/broad/{sample}_peaks.xls",
    log:
        "results/00log/macs2/{sample}_{control}_macs2.log"
    params:
        out_dir       = "results/03peak_macs2/{sample}_{control}/broad",
        macs2_params  = config["params"]["macs2"]["pk_calling"],
        pvalue        = config["params"]["macs2"]["pvalue_broad"],
        pvalue_narrow = config["params"]["macs2"]["pvalue"],
        gsize         = config["params"]["macs2"]["gsize"],
        paired_end    = lambda w: "--format BAM --nomodel" if is_single_end(w.sample) else "--format BAMPE"
    message: 
        "call_peaks macs2 with input for sample {input.case}"
    shell:
        """
        if [[ {wildcards.control} != "NoInput"]]
        then
            macs2 callpeak {params.paired_end} \
                --broad --broad-cutoff {params.pvalue} \
                --treatment {input.case} \
                --control {input.reference} \
                --gsize {params.gsize} \
                --outdir {params.out_dir} \
                --name {wildcards.sample} \
                --pvalue {params.pvalue_narrow} \
                {params.macs2_params} 2> {log}     
        else
            macs2 callpeak {params.paired_end} \
                --broad --broad-cutoff {params.pvalue} \
                --treatment {input.case} \
                --gsize {params.gsize} \
                --outdir {params.out_dir} \
                --name {wildcards.sample} \
                --pvalue {params.pvalue_narrow} \
                {params.macs2_params} 2> {log}
        fi                   
        """


rule filter_peaks:
    input:
        rules.call_peaks.output.narrowPeak
    output:
        bed_filt = "results/03peak_macs2/{sample}_{control}/{sample}_peaks_p{filt}.bed"
    shell:
        """
        awk "\$8 >= {wildcards.filt}" {input} | cut -f1-4,8 > {output.bed_filt}
        """


rule peakAnnot:
    input:
        lambda w: rules.filter_peaks.output.bed_filt if w.pkcaller == "macs2" else rules.seacr_noIgG.output.peaks
    output:
        annot             = "results/04peak_annot/{pkcaller}/{sample}_{control}/{sample}_peaks_p{filt}.annot",
        promo_bed_targets = "results/04peak_annot/{pkcaller}/{sample}_{control}/{sample}_peaks_p{filt}_promoTargets.bed",
        promoTargets      = "results/04peak_annot/{pkcaller}/{sample}_{control}/{sample}_peaks_p{filt}_promoTargets.txt",
        promoBed          = "results/04peak_annot/{pkcaller}/{sample}_{control}/{sample}_peaks_p{filt}_promoPeaks.bed",
        distalBed         = "results/04peak_annot/{pkcaller}/{sample}_{control}/{sample}_peaks_p{filt}_distalPeaks.bed"
    params:
        before = config["promoter"]["bTSS"],
        after  = config["promoter"]["aTSS"],
        genome = lambda wildcards: SAMPLES.GENOME[wildcards.sample]
    log: 
        "results/00log/peakAnnot/{pkcaller}_{sample}_{control}_{filt}_peakanot"
    message:
        "Annotating peaks for {wildcards.sample}"
    shell:
        """
        Rscript --vanilla workflow/scripts/PeakAnnot.R {input} {params.before} {params.after}   \
            {output.annot} {output.promo_bed_targets} {output.promoTargets} {output.promoBed} {output.distalBed} {params.genome}
        """