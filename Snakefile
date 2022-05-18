# If rules are split across multiple snakefiles, list them here
# include: "rules-A"
# include: "rules-B"
from collections import defaultdict

rule all:
    input:
        expand("plots/{which}.{ranges}.{sample}.phaseogram.gene.png", sample=config['reads'], ranges=['126_166', '96_136'], which=['sort','nonrep','rep']),
        # expand('mapping/map.{sample}.ref_asm_mac.nonrep.bam', sample=config['reads']),
        # expand("mapping/map.{sample}.ref_asm_mac.sort.bam", sample=config['reads'])

rule filter_repetitive:
    input:
        bed=config['ref_asm_mac_repbed'],
        bam="mapping/map.{sample}.ref_asm_mac.sort.bam"
    output:
        'mapping/map.{sample}.ref_asm_mac.rep.bam'
    conda: 'envs/mappers.yml'
    log: 'logs/filter_repetitive.{sample}.log'
    threads: 8
    shell:
        r"""
        samtools view -bh -o {output} -L {input.bed} --threads {threads} {input.bam};
        samtools index {output}
        """

rule filter_nonrepetitive:
    input:
        bed='mapping/notrep.bed',
        bam="mapping/map.{sample}.ref_asm_mac.sort.bam"
    output:
        'mapping/map.{sample}.ref_asm_mac.nonrep.bam'
    conda: 'envs/mappers.yml'
    log: 'logs/filter_nonrepetitive.{sample}.log'
    threads: 8
    shell:
        r"""
        samtools view -bh -o {output} -L {input.bed} --threads {threads} {input.bam};
        samtools index {output}
        """

rule nonrepetitive:
    # Get regions that are not annotated as low-complexity repeats
    input:
        bed=config['ref_asm_mac_repbed'],
        g='mapping/ref.ctgsizes'
    output:
        'mapping/notrep.bed'
    conda: 'envs/mappers.yml'
    shell:
        r"""
        bedtools complement -i {input.bed} -g {input.g} > {output}
        """

rule chromsize:
    # contig lengths in genome assembly required for bedtools complement
    input:
        config['ref_asm_mac']
    output:
        'mapping/ref.ctgsizes'
    run:
        ctglens = defaultdict(int)
        with open(input[0], 'r') as fh:
            for line in fh:
                if line.startswith('>'):
                    curr_ctg = line.rstrip()[1:].split()[0] # Split on whitespace
                else:
                    ctglens[curr_ctg] += len(line.rstrip())
        with open(output[0], 'w') as fo:
            for ctg in ctglens:
                fo.write(ctg + '\t' + str(ctglens[ctg]) + '\n')


rule phaseogram:
    input:
        bam="mapping/map.{sample}.ref_asm_mac.{which}.bam",
        gff=config['ref_asm_mac_annot']
    output:
        'plots/{which}.{fragmin}_{fragmax}.{sample}.phaseogram.gene.png'
    log: 'logs/phaseogram.{fragmin}_{fragmax}.{sample}.{which}.log'
    conda: 'mnutils/env.yml'
    params:
        mnutils='workflow/mnutils/mnutils.py',
        prefix='plots/{which}.{fragmin}_{fragmax}.{sample}'
    shell:
        r"""
        python {params.mnutils} -i {input.bam} -o {params.prefix} --min_tlen {wildcards.fragmin} --max_tlen {wildcards.fragmax} --gff {input.gff} --feature gene --phaseogram --dump &> {log}
        """

rule sort_mapping:
    input:
        bam="mapping/map.{sample}.ref_asm_mac.bam",
        ref=config['ref_asm_mac']
    output:
        bam=protected("mapping/map.{sample}.ref_asm_mac.sort.bam"),
        idx="mapping/map.{sample}.ref_asm_mac.sort.bam.bai"
    log: "logs/sort_mapping.{sample}.log"
    threads: 8
    conda: "envs/mappers.yml"
    shell:
        r"""
        samtools sort --threads {threads} --reference {input.ref} -T /tmp/map.{wildcards.sample}.ref_asm_mac.sort -o {output.bam} {input.bam} 2> {log};
        samtools index -@ {threads} {output.bam} 2>> {log};
        """

rule map_pe:
    input:
        fwd=lambda wildcards: config['reads'][wildcards.sample]['fwd'],
        rev=lambda wildcards: config['reads'][wildcards.sample]['rev'],
        ref=config['ref_asm_mac']
    output:
        outm=temp("mapping/map.{sample}.ref_asm_mac.bam")
        # ihist="mapping/map.{sample}.ref_asm_mac.ihist"
    log: "logs/map_pe.{sample}.log"
    conda: "envs/mappers.yml"
    threads: 16
    shell:
        r"""
        minimap2 -ax sr -t {threads} -N 2 {input.ref} {input.fwd} {input.rev} 2> {log} | samtools view -b > {output}
        """
        # "bbmap.sh in={input.fwd} in2={input.rev} threads={threads} ref={input.ref} nodisk semiperfectmode pairlen=1000 pairedonly ihist={output.ihist} outm={output.outm}"

rule merge_reads:
    input:
        fwd=lambda wildcards: config['reads'][wildcards.sample]['fwd'],
        rev=lambda wildcards: config['reads'][wildcards.sample]['rev']
    output:
        merged="reads/{sample}.merged.fastq.gz",
        unmerged="reads/{sample}.unmerged.fastq.gz",
        ihist="reads/{sample}.ihist.txt"
    log: "merge_reads.{sample}.log"
    conda: "envs/mappers.yml"
    threads: 16
    shell:
        "bbmerge.sh threads={threads} in={input.fwd} in2={input.rev} out={output.merged} outu={output.unmerged} ihist={output.ihist} &> {log}"
