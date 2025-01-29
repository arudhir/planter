import json
from pathlib import Path

rule index:
    input:
        transcriptome = rules.rename_headers.output.fasta  # Higher mapping rate because incl. non-coding regions
        # transcriptome = rules.transdecoder.output.longest_orfs_cds  # Lower mapping rate because coding onlyqclear
    params:
        # idx = directory('{outdir}/salmon_index')
        idx = directory(Path(config['outdir']) / '{sample}/salmon_index')
    output:
        done = touch(Path(config['outdir']) / '{sample}/salmon_index/build_index.done')
    log: Path(config['outdir']) / '{sample}/logs/salmon_idx.log'
    threads: workflow.cores * 0.5
    run:
        shell(
            'salmon index '
            '--transcripts {input.transcriptome} '
            '--index {params.idx} '
            '--threads {threads} '
            # '{config[index_params]} '
            '2> {log}'
        )

rule quant:
    input:
        idx = rules.index.output,
        r1 = rules.filter_rrna.output.r1,
        r2 = rules.filter_rrna.output.r2
    output:
        quant_dir = directory(Path(config['outdir']) / '{sample}/quants'),
        tsv = Path(config['outdir']) / '{sample}/quants/quant.sf',
        formatted_tsv = Path(config['outdir']) / '{sample}/quants/{sample}.quant.tsv',
        json = Path(config['outdir']) / '{sample}/quants/{sample}.quant.json',
        stats = Path(config['outdir']) / '{sample}/quants/aux_info/meta_info.json',
#        bam = '{outdir}/quants/{sample}/{sample}.sort.bam',
#        unmapped_names = '{outdir}/quants/{sample}/aux_info/unmapped_names.txt',
#        tmpbam = temp('{outdir}/quants/{sample}/{sample}.bam.tmp'),
    params:
        idx = lambda wildcards: Path(config['outdir']) / f'{wildcards.sample}/salmon_index',
    log:
        Path(config['outdir']) / '{sample}/logs/{sample}_quant.log'
    threads: 8
    run:
        shell(
            'salmon quant '
            '-i {params.idx} '
            '-l A '
            '-1 {input.r1} '
            '-2 {input.r2} '
            '-o {output.quant_dir} '
            '-p {threads} '
            '--writeUnmappedNames '
            '--minAssignedFrags 1 '
            '2> {log}'
        )

        # Add RNAseq_ID to Salmon quants
        import pandas as pd
        df = pd.read_csv(output.tsv, delimiter='\t')

        df['sample'] = wildcards.sample
        df.to_csv(output.formatted_tsv, sep='\t')

        # Make it into a JSON just cuz
        with open(output.json, 'w+') as f:
            json.dump(df.to_dict('records'), f)


