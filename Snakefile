# If rules are split across multiple snakefiles, list them here
# include: "rules-A"
# include: "rules-B"

rule all:
    input:
        expand("reads/{sample}.merged.fastq.gz", sample=config['reads'])

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
