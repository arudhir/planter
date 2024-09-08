rule download_reads:
    output:
        r1 = Path(config['outdir']) / '{sample}/illumina/raw/{sample}_1.fastq',
        r2 = Path(config['outdir']) / '{sample}/illumina/raw/{sample}_2.fastq'
    params:
        outdir = lambda wildcards: Path(config['outdir']) / f'{wildcards.sample}/illumina/raw'
    run:
        shell(
            """
            fastq-dump --split-files {wildcards.sample} -O {params.outdir}
            """
        )

rule compress_reads:
    input:
        r1 = rules.download_reads.output.r1,
        r2 = rules.download_reads.output.r2
    output:
        r1 = Path(config['outdir']) / "{sample}/illumina/raw/{sample}_1.fastq.gz",
        r2 = Path(config['outdir']) / "{sample}/illumina/raw/{sample}_2.fastq.gz"
    run:
        shell(
            """
            pigz {input.r1}
            pigz {input.r2}
            """
        )


rule fastp_raw:
    input:
        r1 = rules.compress_reads.output.r1,
        r2 = rules.compress_reads.output.r2
    output:
        r1 = Path(config['outdir']) / "{sample}/illumina/fastp_raw/{sample}.1.fq.gz",
        r2 = Path(config['outdir']) / "{sample}/illumina/fastp_raw/{sample}.2.fq.gz",
        json = Path(config['outdir']) / "{sample}/illumina/fastp_raw/{sample}_fastp.json",
        html = Path(config['outdir']) / "{sample}/illumina/fastp_raw/{sample}_fastp.html"
    threads: workflow.cores
    run:
        shell(
            'fastp '
            '-i {input.r1} '
            '-I {input.r2} '
            '-o {output.r1} '
            '-O {output.r2} '
            '--json {output.json} '
            '--html {output.html}'
        )


rule filter_rrna:
    input:
        r1 = rules.fastp_raw.output.r1,
        r2 = rules.fastp_raw.output.r2
    output:
        r1 = Path(config['outdir']) / "{sample}/illumina/rrna_filtered/{sample}.1.fq.gz",
        r2 = Path(config['outdir']) / "{sample}/illumina/rrna_filtered/{sample}.2.fq.gz",
        stats = Path(config['outdir']) / "{sample}/illumina/rrna_filtered/{sample}_rRNA_filter.stats"
    log: 
        Path(config['outdir']) / "{sample}/illumina/rrna_filtered/{sample}_rRNA_filter.log"
    threads: workflow.cores * 0.25
    run:
        shell(
            'bbduk.sh '
            'in1={input.r1} '
            'in2={input.r2} '
            'out={output.r1} '
            'out2={output.r2} '
            'refstats={output.stats} '
            'ref={config[rrna_db]} '
            '{config[bbduk_parameters]} '
            '> {log} 2>&1'
        )


# rule trimmomatic:
#     input:
#         r1 = rules.fastp_raw.output.r1,
#         r2 = rules.fastp_raw.output.r2
#     output:
#         r1 = ILMN_READ_DIR / 'trimmed/{sample}_trimmed.1.fq.gz',
#         r1u = ILMN_READ_DIR / 'trimmed/{sample}_trimmed.1u.fq.gz',
#         r2 = ILMN_READ_DIR / 'trimmed/{sample}_trimmed.2.fq.gz',
#         r2u = ILMN_READ_DIR / 'trimmed/{sample}_trimmed.2u.fq.gz'
#     threads: workflow.cores
#     run:
#         shell(
#             'java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar '
#             'PE '
#             '-threads {threads} '
#             '{input.r1} {input.r2} '
#             '{output.r1} {output.r1u} '
#             '{output.r2} {output.r2u} '
#             '{config[trimmomatic_parameters]}'
#         )


rule normalize:
    input:
        r1 = rules.filter_rrna.output.r1,
        r2 = rules.filter_rrna.output.r2
    output:
        r1 = Path(config['outdir']) / "{sample}/illumina/normalized/{sample}_normalized.1.fq.gz",
        r2 = Path(config['outdir']) / "{sample}/illumina/normalized/{sample}_normalized.2.fq.gz"
    params:
        memory = '20g'
    threads: workflow.cores / 2
    run:
        shell(
            'bbnorm.sh '
            'in1={input.r1} '
            'in2={input.r2} '
            'out1={output.r1} '
            'out2={output.r2} '
            'target=200 '
            'min=5 '
            'threads={threads} '
            '-Xmx{params.memory}'
        )


rule fastp_processed:
    input:
        r1 = rules.normalize.output.r1,
        r2 = rules.normalize.output.r2
    output:
        r1 = (Path(config['outdir']) / "{sample}/illumina/processed/{sample}.1.fq.gz"),
        r2 = (Path(config['outdir']) / "{sample}/illumina/processed/{sample}.2.fq.gz"),
        json = (Path(config['outdir']) / "{sample}/illumina/processed/{sample}_fastp.json"),
        html = (Path(config['outdir']) / "{sample}/illumina/processed/{sample}_fastp.html")
    threads: workflow.cores
    run:
        shell(
            'fastp '
            '-i {input.r1} '
            '-I {input.r2} '
            '-o {output.r1} '
            '-O {output.r2} '
            '--json {output.json} '
            '--html {output.html}'
        )
