
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
        # bam = lambda wc: expand("bam/" + wc.sample + "/selected/{bam}.bam", bam = bam_per_sample[wc.sample]) if wc.sample in bam_per_sample else "FOOBAR",
        # bai = lambda wc: expand("bam/" + wc.sample + "/selected/{bam}.bam.bai", bam = bam_per_sample[wc.sample]) if wc.sample in bam_per_sample else "FOOBAR",
        config['input_bam_location'] + "{sample}/selected/{bam}.bam"
        # excl = "log/exclude_file"
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
            -w {wildcards.window} \
            {input.bam} 
        """

