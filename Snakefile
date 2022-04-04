# If rules are split across multiple snakefiles, list them here
# include: "rules-A"
# include: "rules-B"

rule all:
    input:
        expand("mapping/map.{sample}.ref_asm_mac.bam", sample=config['reads'])

rule map_pe:
    input:
        fwd=lambda wildcards: config['reads'][wildcards.sample]['fwd'],
        rev=lambda wildcards: config['reads'][wildcards.sample]['rev'],
        ref=config['ref_asm_mac']
    output:
        outm="mapping/map.{sample}.ref_asm_mac.bam"
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
