################################################################################
# Base
# Base
FROM phusion/baseimage:noble-1.0.0
ENV TOOLS=/tools
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    autotools-dev \
    automake \
    awscli \
    bamtools \
    bioperl \ 
    cmake \
    coinor-cbc \
    curl \
    git \
    gnuplot \
    graphviz \
    libbamtools-dev \
    libboost-all-dev \
    libdatetime-perl \
    libdigest-md5-perl \
    liblzma-dev \
    libncurses5-dev \
    libncursesw5-dev \
    libtbbmalloc2 \
    libxml-simple-perl \
    libbz2-dev \
    nodejs \
    openjdk-11-jre-headless \
    parallel \
    pigz \
    python3 \
    python3-pip \
    python-is-python3 \
    python3-venv \
    pipx \
    progressivemauve \
    fastqc \
    wget \
    salmon

# Install Rust
RUN curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs -y | sh

# Install uv
RUN pip install uv --break-system-packages

# TODO: Move to uv Install uv
#RUN curl -LsSf https://astral.sh/uv/install.sh | sh
# RUN pip install uv --break-system-packages
# RUN uv init
# RUN uv venv
# RUN source .venv/bin/activate
# #COPY --from=ghcr.io/astral-sh/uv:0.4.6 /uv /bin/uv

# # Set up Python virtual environment
# ENV VIRTUAL_ENV=/opt/venv
# RUN python3 -m venv $VIRTUAL_ENV
# ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Verify we're using the correct Python and pip
# RUN which python3 && which pip

# Copy requirements file
# COPY requirements.txt .

# Install Python packages
# RUN pip install --upgrade pip && pip install -r requirements.txt

# Set the working directory
WORKDIR $TOOLS

#################################################################################
### 3rd party tools
#################################################################################

###########
## BBMap ##
###########
RUN wget -O BBMap_37.41.tar.gz https://sourceforge.net/projects/bbmap/files/BBMap_37.41.tar.gz/download \
    && gunzip BBMap_37.41.tar.gz \
    && tar -xvf BBMap_37.41.tar

#####################
## Samtools et. al ##
#####################
RUN wget https://github.com/samtools/samtools/releases/download/1.9/samtools-1.9.tar.bz2 \
    && tar xvf samtools-1.9.tar.bz2 \
    && cd samtools-1.9 \
    && ./configure \
    && make \
    && make install
ENV SAMTOOLS="${TOOLS}/samtools-1.9/samtools"

RUN wget https://github.com/samtools/bcftools/releases/download/1.9/bcftools-1.9.tar.bz2 \
    && tar xvf bcftools-1.9.tar.bz2 \
    && cd bcftools-1.9 \
    && ./configure \
    && make \
    && make install
ENV BCFTOOLS="${TOOLS}/bcftools-1.9/bcftools"

RUN wget https://github.com/samtools/htslib/releases/download/1.9/htslib-1.9.tar.bz2 \
    && tar xvf htslib-1.9.tar.bz2 \
    && cd htslib-1.9 \
    && ./configure \
    && make \
    && make install
ENV HTSLIB="${TOOLS}/htslib-1.9/htslib"

############
### BLAST ##
############
RUN wget ftp://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/2.9.0/ncbi-blast-2.9.0+-x64-linux.tar.gz \
    && tar -zxvf ncbi-blast-2.9.0+-x64-linux.tar.gz

##############
### tbl2asn ##
##############
## More info here: https://www.ncbi.nlm.nih.gov/genbank/tbl2asn2/
# Update: tbl2asn no longer exists, replaced by table2asn: https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/
RUN wget https://ftp.ncbi.nlm.nih.gov/asn1-converters/by_program/table2asn/linux64.table2asn.gz \
    && gunzip linux64.table2asn.gz \
    && mv linux64.table2asn table2asn \
    && chmod a+x table2asn

##################
###### seqtk #####
##################
#RUN git clone https://github.com/lh3/seqtk.git \
    #&& cd seqtk \
    #&& make

##################
###### fastp #####
##################
RUN wget http://opengene.org/fastp/fastp \
    && chmod a+x ./fastp

###################
###### SPAdes #####
###################
WORKDIR ${TOOLS}
RUN wget https://github.com/ablab/spades/releases/download/v4.0.0/SPAdes-4.0.0-Linux.tar.gz \
    && tar -xvzf "${TOOLS}/SPAdes-4.0.0-Linux.tar.gz" \
    && rm "${TOOLS}/SPAdes-4.0.0-Linux.tar.gz" 
    # && cp "{TOOLS}/SPAdes-3.14-Linux.tar.gz" 

################
#### MMseqs2 ###
################
#WORKDIR ${TOOLS}
#RUN wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz \
    #&& tar xvfz mmseqs-linux-avx2.tar.gz 

###############
#### csvtk ####
###############
#WORKDIR ${TOOLS}
#RUN wget https://github.com/shenwei356/csvtk/releases/download/v0.23.0/csvtk_linux_amd64.tar.gz \
    #&& tar xvf csvtk_linux_amd64.tar.gz

####################
### SRA Toolkit ####
####################
WORKDIR ${TOOLS}
RUN wget https://ftp-trace.ncbi.nlm.nih.gov/sra/sdk/3.1.1/sratoolkit.3.1.1-ubuntu64.tar.gz \
    && tar -zxvf sratoolkit.3.1.1-ubuntu64.tar.gz

######################
#### TransDecoder ####
######################
WORKDIR ${TOOLS}
RUN wget https://github.com/TransDecoder/TransDecoder/archive/refs/tags/TransDecoder-v5.7.1.tar.gz \
    && tar -zxvf TransDecoder-v5.7.1.tar.gz

################
#### eggNOG ####
################
WORKDIR ${TOOLS}
RUN git clone https://github.com/eggnogdb/eggnog-mapper

################
#### mmseqs ####
################
WORKDIR ${TOOLS}
RUN wget https://mmseqs.com/latest/mmseqs-linux-avx2.tar.gz \
    && tar xvfz mmseqs-linux-avx2.tar.gz

################
#### seqkit ####
################
WORKDIR ${TOOLS}
RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.9.0/seqkit_linux_amd64.tar.gz \
    && tar xvf seqkit_linux_amd64.tar.gz

###############
### DuckDB ###
###############
WORKDIR ${TOOLS}
RUN wget https://github.com/duckdb/duckdb/releases/download/v1.2.0/duckdb_cli-linux-amd64.zip \
    && unzip duckdb_cli-linux-amd64.zip \
    && rm duckdb_cli-linux-amd64.zip

### Set PATH and do some extraneous steps
ENV PATH=/:/usr/src/planter/planter:/tools/ntHits-ntHits-v0.0.1/:/tools/minimap2:/tools/racon/build/bin:/tools/ntEdit/:/tools/racon-v1.3.1/build/bin:/tools/augustus-3.3.2/bin:/tools/augustus-3.3.2/bin/scripts:/tools/ncbi-blast-2.9.0+/bin:/tools/Flye/bin:/tools/prokka-1.14.0/bin/:/tools/barrnap-0.8/bin:/tools/bbmap:/tools/bowtie2-2.3.2/:/tools/SPAdes-4.0.0-Linux/bin:/tools/Filtlong/bin:/tools/canu-1.8/Linux-amd64/bin:/tools:/usr/bin:/tools/SKESA:/tools/miniasm/:/tools/seqtk:/tools/quast-5.0.2:/tools/minigraph:/tools/SPAdes-3.13.0-Linux/bin:$PATH:/tools/mmseqs/bin/:/tools/sratoolkit.3.1.1-ubuntu64/bin:/tools/TransDecoder-TransDecoder-v5.7.1/:/tools/eggnog-mapper:/tools/mmseqs/bin/:$PATH

#RUN rm ${TOOLS}/*.gz && \
    #rm ${TOOLS}/*.zip

#################################################################################
# Workflow and homemade scripts
# Python environment setup with uv
ENV APP_HOME=/usr/src/planter
WORKDIR $APP_HOME

# Install uv properly
RUN curl -LsSf https://astral.sh/uv/install.sh | sh

# Create virtual environment
ENV VIRTUAL_ENV=/opt/venv
RUN uv venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Copy Python package files
COPY pyproject.toml .
COPY planter/ planter/
# COPY workflow/ workflow/
COPY requirements.txt .

# Install dependencies and package
RUN uv pip install -r requirements.txt
RUN uv pip install -e .

# Patch Snakemake workflow.py to fix logger issues
RUN python -c "import os, re; \
workflow_path = '/opt/venv/lib/python3.12/site-packages/snakemake/workflow.py'; \
content = open(workflow_path, 'r').read(); \
patched = content \
.replace('logger.run_info(', 'logger.info(\"Run info: \", ') \
.replace('logger.resources_info(', 'logger.info(\"Resources info: \", ') \
.replace('logger.host_info(', 'logger.info(\"Host info: \", ') \
.replace('logger.logfile_hint(', 'logger.info(\"Logfile hint: \", ') \
.replace('logger.get_logfile()', 'None'); \
# Find and replace logger.info() calls that don't have a message parameter
patched = re.sub(r'logger\.info\(\s*\)', 'logger.info(\"\")', patched); \
# Comment out any remaining logfile_hint calls
patched = re.sub(r'logger\.logfile_hint\(.*?\)', '# logger.logfile_hint call removed', patched); \
open(workflow_path, 'w').write(patched); \
print('Patched Snakemake workflow.py successfully')"

# Create necessary directories
RUN mkdir -p $APP_HOME/inputs $APP_HOME/outputs

COPY .snakemake .snakemake
