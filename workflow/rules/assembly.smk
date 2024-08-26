TRINITY_OUT_DIR = Path(config['outdir']) / 'trinity'
RNASPADES_OUT_DIR = Path(config['outdir']) / 'rnaspades'
TRANSDECODER_OUT_DIR = Path(config['outdir']) / 'transdecoder'

rule trinity:
    input:
        r1 = rules.fastp_processed.output.r1,
        r2 = rules.fastp_processed.output.r2
    output:
        trinity_outdir = directory(TRINITY_OUT_DIR / '{sample}_trinity'),
        fasta = TRINITY_OUT_DIR / '{sample}_trinity/Trinity.tmp.fasta',
    threads: workflow.cores
    run:
        shell(
            'Trinity '
            '--seqType fq '
            '--max_memory 25G '
            '--output {output.trinity_outdir} '
            '--left {input.r1} '
            '--right {input.r2} '
        )

rule rnaspades:
    input:
        r1 = rules.fastp_processed.output.r1,
        r2 = rules.fastp_processed.output.r2
    output:
        rnaspades_outdir = directory(RNASPADES_OUT_DIR / '{sample}_rnaspades'),
        fasta = RNASPADES_OUT_DIR / '{sample}_rnaspades/{sample}_transcripts.fasta'
    params:
        original_fasta_filename = RNASPADES_OUT_DIR / '{sample}/transcripts.fasta',
        memory = 25
    threads: workflow.cores
    run:
        shell(
            'rnaspades.py '
            '-1 {input.r1} '
            '-2 {input.r2} '
            '-o {output.rnaspades_outdir} '
            '--threads {threads} '
            '--memory {params.memory}'
        )

        shell('mv {params.original_fasta_filename} {output.fasta}')


rule transdecoder:
    input:
        fasta = rules.rnaspades.output.fasta
    output:
        outdir = directory(TRANSDECODER_OUT_DIR / '{sample}_transdecoder'),
        longest_orfs = (TRANSDECODER_OUT_DIR / '{sample}_transdecoder') / 'longest_orfs.cds',
        longest_orfs_gff3 = (TRANSDECODER_OUT_DIR / '{sample}_transdecoder') / 'longest_orfs.gff3',
        longest_orfs_pep = (TRANSDECODER_OUT_DIR / '{sample}_transdecoder') / 'longest_orfs.pep',
    run:
        shell(
            """
            TransDecoder.LongOrfs -t {input.fasta} -O {output.outdir}
            """
        )

        shell(
            """
            TransDecoder.Predict -t {input.fasta} --no_refine_starts -O {output.outdir}
            """
        )


#rule busco:
#    pass
#
