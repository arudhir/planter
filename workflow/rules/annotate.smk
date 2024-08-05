EGGNOG_OUT_DIR = Path(config['outdir']) / 'eggnog'

rule eggnog:
    input:
        transcriptome = rules.transdecoder.output.longest_orfs_pep
    output:
        annotations = EGGNOG_OUT_DIR / '{sample}_eggnog.emapper.annotations.xlsx'
    threads: 16
    params:
        eggnog_datadir = '/mnt/data',
        eggnog_outdir = directory(EGGNOG_OUT_DIR / '{sample}_eggnog'),
    run:
        shell(
            '/opt/eggnog-mapper/emapper.py '
            '-i {input.transcriptome} '
            '--itype proteins '
            '--output {params.eggnog_outdir} '
            '--cpu {threads} '
            '--data_dir {params.eggnog_datadir} '
            '--excel '
            # '--output_dir {output.eggnog_outdir} '
            '--override'
        )
