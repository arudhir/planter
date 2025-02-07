import shutil
from pathlib import Path

TRINITY_OUT_DIR = Path(config['outdir']) / 'trinity'
RNASPADES_OUT_DIR = Path(config['outdir']) / 'rnaspades'
TRANSDECODER_OUT_DIR = Path(config['outdir']) / 'transdecoder'

rule trinity:
    input:
        r1 = rules.fastp_processed.output.r1,
        r2 = rules.fastp_processed.output.r2
    output:
        trinity_outdir = directory(Path(config['outdir']) / '{sample}/trinity'),
        fasta = Path(config['outdir']) / '{sample}/trinity/Trinity.tmp.fasta',
    threads: workflow.cores
    run:
        shell(
            'Trinity '
            '--seqType fq '
            '--max_memory 30G '
            '--output {output.trinity_outdir} '
            '--left {input.r1} '
            '--right {input.r2} '
        )

rule rnaspades:
    input:
        r1 = rules.fastp_processed.output.r1,
        r2 = rules.fastp_processed.output.r2
    output:
        fasta = Path(config['outdir']) / '{sample}/rnaspades/transcripts.fasta'
    params:
        memory = 20,
        rnaspades_outdir = directory(Path(config['outdir']) / '{sample}/rnaspades'),
    threads: workflow.cores / 2
    run:
        shell(
            'rnaspades.py '
            '-1 {input.r1} '
            '-2 {input.r2} '
            '-o {params.rnaspades_outdir} '
            '--checkpoints last '
            '--threads {threads} '
            '--memory {params.memory}'
        )

rule rename_assembly:
    input:
        original = rules.rnaspades.output.fasta
    output:
        renamed = Path(config['outdir']) / '{sample}/rnaspades/{sample}_transcripts.fasta'
    run:
        shutil.copy(input.original, output.renamed)

rule rename_headers:
    input:
        fasta = rules.rename_assembly.output.renamed
    output:
        fasta = Path(config['outdir']) / '{sample}/rnaspades/{sample}_transcripts_renamed.fasta'
    shell:
        "./planter/scripts/seqhash_rename.py --input {input.fasta} --output {output.fasta}"

rule transdecoder:
    input:
        fasta = rules.rename_headers.output.fasta
    output:
        outdir = directory(Path(config['outdir']) / '{sample}/transdecoder'),
        longest_orfs_cds = (Path(config['outdir']) / '{sample}/transdecoder') / '{sample}.cds',
        longest_orfs_gff3 = (Path(config['outdir']) / '{sample}/transdecoder') / '{sample}.gff3',
        longest_orfs_pep = (Path(config['outdir']) / '{sample}/transdecoder') / '{sample}.pep',
        longest_orfs_bed = (Path(config['outdir']) / '{sample}/transdecoder') / '{sample}.bed',
    params:
        original_cds = lambda wildcards: Path(config['outdir']) / '{wildcards.sample}/transdecoder/{wildcards.sample}_transcripts_renamed.fasta.transdecoder.cds',
        original_gff3 = lambda wildcards: Path(config['outdir']) / '{wildcards.sample}/transdecoder/{wildcards.sample}_transcripts_renamed.fasta.transdecoder.gff3',
        original_pep = lambda wildcards: Path(config['outdir']) / '{wildcards.sample}/transdecoder/{wildcards.sample}_transcripts_renamed.fasta.transdecoder.pep',
        original_bed = lambda wildcards: Path(config['outdir']) / '{wildcards.sample}/transdecoder/{wildcards.sample}_transcripts_renamed.fasta.transdecoder.bed'
    run:
        shell(
            """
            TransDecoder.LongOrfs -t {input.fasta} --output_dir {output.outdir} --complete_orfs_only
            """
        )

        shell(
            """
            TransDecoder.Predict -t {input.fasta} --no_refine_starts --output_dir {output.outdir}
            """
        )
        shell(
            f"""
            mv {params.original_bed} {output.longest_orfs_bed}
            mv {params.original_cds} {output.longest_orfs_cds}
            mv {params.original_gff3} {output.longest_orfs_gff3}
            mv {params.original_pep} {output.longest_orfs_pep}
            """
        )


#rule busco:
#    pass
#
