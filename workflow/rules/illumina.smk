ILMN_READ_DIR = Path(config['outdir']) / 'illumina'

rule download_reads:
    output:
        r1 = ILMN_READ_DIR / "raw/{sample}_1.fastq",
        r2 = ILMN_READ_DIR / "raw/{sample}_2.fastq"
    params:
        outdir = ILMN_READ_DIR / 'raw'
    run:
        shell(
            """
            fasterq-dump --split-files {wildcards.sample} -O {params.outdir}
            """
        )

rule compress_reads:
    input:
        r1 = rules.download_reads.output.r1,
        r2 = rules.download_reads.output.r2
    output:
        r1 = ILMN_READ_DIR / "raw/{sample}_1.fastq.gz",
        r2 = ILMN_READ_DIR / "raw/{sample}_2.fastq.gz"
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
        r1 = (ILMN_READ_DIR / 'fastp_raw/{sample}.1.fq.gz'),
        r2 = (ILMN_READ_DIR / 'fastp_raw/{sample}.2.fq.gz'),
        json = (ILMN_READ_DIR / 'fastp_raw/{sample}_fastp.json'),
        html = (ILMN_READ_DIR / 'fastp_raw/{sample}_fastp.html')
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

rule trimmomatic:
    input:
        r1 = rules.fastp_raw.output.r1,
        r2 = rules.fastp_raw.output.r2
    output:
        r1 = ILMN_READ_DIR / 'trimmed/{sample}_trimmed.1.fq.gz',
        r1u = ILMN_READ_DIR / 'trimmed/{sample}_trimmed.1u.fq.gz',
        r2 = ILMN_READ_DIR / 'trimmed/{sample}_trimmed.2.fq.gz',
        r2u = ILMN_READ_DIR / 'trimmed/{sample}_trimmed.2u.fq.gz'
    threads: 8
    run:
        shell(
            'java -jar /opt/Trimmomatic-0.39/trimmomatic-0.39.jar '
            'PE '
            '-threads {threads} '
            '{input.r1} {input.r2} '
            '{output.r1} {output.r1u} '
            '{output.r2} {output.r2u} '
            '{config[trimmomatic_parameters]}'
        )

rule fastp_trimmed:
    input:
        r1 = rules.trimmomatic.output.r1,
        r2 = rules.trimmomatic.output.r2
    output:
        r1 = (ILMN_READ_DIR / 'fastp_trimmed/{sample}.1.fq.gz'),
        r2 = (ILMN_READ_DIR / 'fastp_trimmed/{sample}.2.fq.gz'),
        json = (ILMN_READ_DIR / 'fastp_trimmed/{sample}_fastp.json'),
        html = (ILMN_READ_DIR / 'fastp_trimmed/{sample}_fastp.html')
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