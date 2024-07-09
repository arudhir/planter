TRINITY_OUT_DIR = Path(config['outdir']) / 'trinity'

rule trinity:
    input:
        r1 = rules.trimmomatic.output.r1,
        r2 = rules.trimmomatic.output.r2
    output:
        trinity_outdir = directory(TRINITY_OUT_DIR / '{sample}_trinity'),
        trinity_fasta = TRINITY_OUT_DIR / '{sample}_trinity/Trinity.tmp.fasta'
    run:
        shell(
            'Trinity '
            '--seqType fq '
            '--max_memory 50G '
            '--output {output.trinity_outdir} '
            '--left {input.r1} '
            '--right {input.r2} '
        )

#rule transdecoder:
#    pass
#
#rule busco:
#    pass
#
