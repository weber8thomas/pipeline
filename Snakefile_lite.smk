# Python imports
import math
from collections import defaultdict
import os, sys

# Snakemake config file
configfile: "Snake.config.yaml"



"""
From Snakemake config file, generate input files architecture : 1 sample per folder and 1 or multiple BAM files per sample
"""

# Retrieve samples and bam files into python lists
sample_list, bam_list = glob_wildcards(config['input_bam_location'] + "{sample}/selected/{bam}.bam")

print(sample_list, type(sample_list))
print(bam_list, type(bam_list))

# cell_per_sample= defaultdict(list)
bam_per_sample = defaultdict(list)

for sample,bam in zip(sample_list,bam_list):
    bam_per_sample[sample].append(bam)
    # cell_per_sample[sample].append(bam.replace('.sort.mdup',''))

# allbams_per_sample = defaultdict(list)
# for sample in samples:
#     allbams_per_sample[sample] = glob_wildcards(config['input_bam_location'] + "{}/all/{{bam}}.bam".format(sample)).bam

# print("Detected {} samples:".format(len(samples)))
# for s in samples:
#     print("  {}:\t{} cells\t {} selected cells".format(s, len(allbams_per_sample[s]), len(bam_per_sample[s])))


import os.path

# Current state of the pipeline:
# ==============================
# * count reads in the bam files (in fixed and variable-width bins of various sizes)
# * determine strand states of each chromosome in each single cell, including SCEs
# * plot all single cell libraries in different window sizes
# * calculate a segmentation into potential SVs using Mosaicatcher



rule all:
    input:
        # expand("counts/{sample}/{window}_fixed.txt.gz", sample=samples, window=[100000])
        "log/exclude_file.temp" 

################################################################################
# Read counting                                                                #
################################################################################

rule generate_exclude_file_1:
    output:
        temp("log/exclude_file.temp")
    input:
        bam = expand(config['input_bam_location'] + "{sample}/selected/{bam}.bam", sample = sample_list[0], bam = bam_per_sample[sample_list[0]][0])
    log:
        "log/generate_exclude_file_1.log"
    params:
        samtools = config["samtools"]
    shell:
        """
        {params.samtools} view -H {input.bam} | awk "/^@SQ/" > {output} 
        """

rule generate_exclude_file_2:
    output:
        "log/exclude_file"
    input:
        "log/exclude_file.temp"
    params:
        chroms = config["chromosomes"]
    run:
        with open(input[0]) as f:
            with open(output[0],"w") as out:
                for line in f:
                    contig = line.strip().split()[1]
                    contig = contig[3:]
                    # if contig not in params.chroms:
                        # print(contig, file = out)


rule mosaic_count_fixed:
    input:
        bam = lambda wc: expand("bam/" + wc.sample + "/selected/{bam}.bam", bam = bam_per_sample[wc.sample]) if wc.sample in bam_per_sample else "FOOBAR",
        bai = lambda wc: expand("bam/" + wc.sample + "/selected/{bam}.bam.bai", bam = bam_per_sample[wc.sample]) if wc.sample in bam_per_sample else "FOOBAR",
        excl = "log/exclude_file"
    output:
        counts = "counts/{sample}/{window}_fixed.txt.gz",
        info   = "counts/{sample}/{window}_fixed.info"
    log:
        "log/{sample}/mosaic_count_fixed.{window}.log"
    params:
        mc_command = config["mosaicatcher"]
    shell:
        """
        echo mosaic_count_fixed && 
        {params.mc_command} count \
            --verbose \
            --do-not-blacklist-hmm \
            -o {output.counts} \
            -i {output.info} \
            -x {input.excl} \
            -w {wildcards.window} \
            {input.bam} 
        """

# rule extract_single_cell_counts:
#     input:
#         "counts/{sample}/{window}_{file_name}.txt.gz"
#     output:
#         "counts-per-cell/{sample}/{cell}/{window,[0-9]+}_{file_name}.txt.gz"
#     shell:
#         """
#         # Issue #1022 (https://bitbucket.org/snakemake/snakemake/issues/1022)
#         zcat {input} | awk -v name={wildcards.cell} -f utils/command1.awk | gzip > {output}
#         """


# ################################################################################
# # Plots                                                                        #
# ################################################################################

# rule plot_mosaic_counts:
#     input:
#         counts = "counts/{sample}/{file_name}.txt.gz",
#         info   = "counts/{sample}/{file_name}.info"
#     output:
#         "plots/{sample}/{file_name}.pdf"
#     log:
#         "log/plot_mosaic_counts/{sample}/{file_name}.log"
#     params:
#         plot_command = "Rscript " + config["plot_script"]
#     shell:
#         """
#         {params.plot_command} {input.counts} {input.info} {outp
#         """

