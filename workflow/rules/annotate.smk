EGGNOG_OUT_DIR = Path(config['outdir']) / 'eggnog'

rule eggnog:
    input:
        transcriptome = rules.transdecoder.output.longest_orfs_pep
    output:
        annotations = Path(config['outdir']) / '{sample}/eggnog/{sample}.emapper.annotations'
    threads: 16
    params:
        eggnog_datadir = '/mnt/data',
        eggnog_outdir = directory(Path(config['outdir']) / '{sample}/eggnog'),
    run:
        shell(
            'source $VIRTUAL_ENV/bin/activate && /tools/eggnog-mapper/emapper.py '
            '-i {input.transcriptome} '
            '--itype proteins '
            '--output {wildcards.sample} '
            '--cpu {threads} '
            '--data_dir {params.eggnog_datadir} '
            '--excel '
            '--output_dir {params.eggnog_outdir} '
            '--report_orthologs '
            '--report_no_hits '
            '--override'
        )

rule analyze_eggnog:
    input:
        annotations = rules.eggnog.output.annotations
    output:
        directory(Path(config['outdir']) / '{sample}/eggnog/plots')
    run:
        shell(
            './planter/scripts/parse_eggnog.py {input.annotations} {output}'
        )
