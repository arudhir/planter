TRINITY_OUT_DIR = Path(config['outdir']) / 'trinity'
TRANSDECODER_OUT_DIR = Path(config['outdir']) / 'transdecoder'

rule trinity:
    input:
        r1 = rules.trimmomatic.output.r1,
        r2 = rules.trimmomatic.output.r2
    output:
        trinity_outdir = directory(TRINITY_OUT_DIR / '{sample}_trinity'),
        cds = TRINITY_OUT_DIR / '{sample}_trinity/Trinity.tmp.fasta',
    threads: 8
    run:
        shell(
            'Trinity '
            '--seqType fq '
            '--max_memory 50G '
            '--output {output.trinity_outdir} '
            '--left {input.r1} '
            '--right {input.r2} '
        )

rule transdecoder:
    input:
        trinity_fasta = rules.trinity.output.cds
    output:
        outdir = directory(TRANSDECODER_OUT_DIR / '{sample}_transdecoder'),
        longest_orfs = (TRANSDECODER_OUT_DIR / '{sample}_transdecoder') / 'longest_orfs.cds',
        longest_orfs_gff3 = (TRANSDECODER_OUT_DIR / '{sample}_transdecoder') / 'longest_orfs.gff3',
        longest_orfs_pep = (TRANSDECODER_OUT_DIR / '{sample}_transdecoder') / 'longest_orfs.pep',
    run:
        shell(
            """
            TransDecoder.LongOrfs -t {input.trinity_fasta} -O {output.outdir}
            """
        )

        shell(
            """
            TransDecoder.Predict -t {input.trinity_fasta} --no_refine_starts -O {output.outdir}
            """
        )


#rule busco:
#    pass
#
