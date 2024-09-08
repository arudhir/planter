################################################################################
# Base
# Base
FROM phusion/baseimage:noble-1.0.0
ENV TOOLS=/tools
ENV DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y \
    autotools-dev \
    automake \
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
    progressivemauve \
    fastqc \
    wget

# TODO: Move to uv Install uv
#RUN curl -LsSf https://astral.sh/uv/install.sh | sh
RUN pip install uv --break-system-packages
#COPY --from=ghcr.io/astral-sh/uv:0.4.6 /uv /bin/uv

# Set up Python virtual environment
ENV VIRTUAL_ENV=/opt/venv
RUN python3 -m venv $VIRTUAL_ENV
ENV PATH="$VIRTUAL_ENV/bin:$PATH"

# Verify we're using the correct Python and pip
RUN which python3 && which pip

# Copy requirements file
COPY requirements.txt .

# Install Python packages
RUN pip install --upgrade pip && pip install -r requirements.txt

# Set the working directory
WORKDIR $TOOLS

#################################################################################
### 3rd party tools
#################################################################################
#ENV TOOLS=/tools
#WORKDIR ${TOOLS}

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

##############
### QUAST ####
##############
#RUN wget https://downloads.sourceforge.net/project/quast/quast-5.0.2.tar.gz \
    #&& tar -xzf quast-5.0.2.tar.gz

##############
### Prokka ###
##############
RUN wget -O prokka-v1.14.tar.gz https://github.com/tseemann/prokka/archive/v1.14.0.tar.gz \
    && tar -xzvf prokka-v1.14.tar.gz

##############
#### Canu ####
##############
RUN wget https://github.com/marbl/canu/releases/download/v1.8/canu-1.8.Linux-amd64.tar.xz \
    && xz -dc canu-1.8.*.tar.xz | tar -xf -

##############
#### Flye ####
##############
RUN git clone https://github.com/fenderglass/Flye \
    && cd Flye \
    && python3 setup.py build \
    && python3 setup.py install

#################
#### Filtlong ###
#################
#RUN git clone https://github.com/rrwick/Filtlong.git \
    #&& cd Filtlong \
    #&& make -j

##################
#### Unicycler ###
##################
#RUN git clone https://github.com/rrwick/Unicycler.git \
    #&& cd Unicycler \
    #&& python3 setup.py install

##################
#### Barrnap #####
##################
#RUN wget -O barrnap-0.8.tar.gz https://github.com/tseemann/barrnap/archive/0.8.tar.gz \
    #&& tar -zxvf barrnap-0.8.tar.gz

##################
#### Bandage #####
##################
##RUN wget https://github.com/rrwick/Bandage/releases/download/v0.8.1/Bandage_Ubuntu_static_v0_8_1.zip \
    ##&& unzip Bandage_Ubuntu_static_v0_8_1.zip

##################
##### seqkit #####
##################
#RUN wget https://github.com/shenwei356/seqkit/releases/download/v2.1.0/seqkit_linux_amd64.tar.gz \
    #&& tar -zxvf seqkit_linux_amd64.tar.gz

##################
##### Pilon ######
##################
#RUN wget https://github.com/broadinstitute/pilon/releases/download/v1.22/pilon-1.22.jar

##################
##### Bowtie2 ####
##################
RUN wget https://sourceforge.net/projects/bowtie-bio/files/bowtie2/2.3.2/bowtie2-2.3.2-linux-x86_64.zip \
    && unzip bowtie2-2.3.2-linux-x86_64.zip

##################
##### Porechop ###
##################
#RUN git clone https://github.com/rrwick/Porechop.git \
    #&& cd Porechop \
    #&& python3 setup.py install

##################
##### Prodigal ###
##################
#RUN wget --output-document prodigal.tar.gz https://github.com/hyattpd/Prodigal/archive/v2.6.3.tar.gz \
    #&& tar -zxvf prodigal.tar.gz \
    #&& cd Prodigal-2.6.3/ \
    #&& make install INSTALLDIR=/usr/bin/

##################
##### Sepp #######
##################
#RUN git clone  https://github.com/smirarab/sepp.git \
    #&& cd sepp \
    #&& python3 setup.py config \
    #&& python3 setup.py install

##################
##### HMMER ######
##################
#RUN wget http://eddylab.org/software/hmmer3/3.1b2/hmmer-3.1b2-linux-intel-x86_64.tar.gz \
    #&& gunzip hmmer-3.1b2-linux-intel-x86_64.tar.gz \
    #&& tar -xf hmmer-3.1b2-linux-intel-x86_64.tar

##################
##### Racon ######
##################
#RUN wget https://github.com/isovic/racon/releases/download/1.3.1/racon-v1.3.1.tar.gz \
    #&& tar -xf racon-v1.3.1.tar.gz \
    #&& cd racon-v1.3.1 \
    #&& mkdir build \
    #&& cd build \
    #&& cmake -DCMAKE_BUILD_TYPE=Release .. \
    #&& make

##################
##### ntEdit #####
##################
#RUN git clone https://github.com/bcgsc/ntEdit.git \
    #&& cd ntEdit \
    #&& make ntedit \
    #&& chmod +x ntedit

##################
##### ntHits #####
##################
#RUN wget https://github.com/bcgsc/ntHits/archive/ntHits-v0.0.1.tar.gz \
    #&& tar -xzvf ntHits-v0.0.1.tar.gz \
    #&& cd ntHits-ntHits-v0.0.1 && ./autogen.sh \
    #&& ./configure && make && make install

##################
##### Minimap2 ###
##################
#WORKDIR $TOOLS
#RUN git clone https://github.com/lh3/minimap2 \
    #&& cd minimap2 \
    #&& make

##################
#### minigraph ###
##################
RUN git clone https://github.com/lh3/minigraph \
    && cd minigraph && make

##################
##### Miniasm ####
##################
WORKDIR $TOOLS
RUN git clone https://github.com/lh3/miniasm  \
    && cd miniasm  \
    && make

##################
##### MUMmer #####
##################
#RUN wget https://github.com/mummer4/mummer/releases/download/v4.0.0beta2/mummer-4.0.0beta2.tar.gz \
    #&& tar xvf mummer-4.0.0beta2.tar.gz \
    #&& cd mummer-4.0.0beta2 \
    #&& ./configure prefix=/usr/bin \
    #&& make \
    #&& make install

##################
##### Shasta #####
##################
#RUN curl -O -L https://github.com/chanzuckerberg/shasta/releases/download/0.6.0/shasta-Linux-0.6.0 \
    #&& chmod ugo+x shasta-Linux-0.6.0

################
#### SKESA #####
################
WORKDIR ${TOOLS}
RUN git clone https://github.com/ncbi/SKESA \
    && cd SKESA \
    && make -f Makefile.nongs

################
#### DIAMOND ###
################
#RUN wget https://github.com/bbuchfink/diamond/releases/download/v2.0.5/diamond-linux64.tar.gz && \
    #tar -xvf diamond-linux64.tar.gz && \
    #rm diamond-linux64.tar.gz

##################
#### Botowraps ###
##################
##RUN git clone git@git.ginkgobioworks.com:ngs-analysis/botowraps.git \
    ##&& cd botowraps \
    ##&& pip3 install .

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

###################
###### pauvre #####
###################
##WORKDIR ${TOOLS}
##RUN git clone https://github.com/conchoecia/pauvre.git \
    ##&& cd ./pauvre \
    ##&& pip3 install biopython==1.74 .

##################
###### Busco #####
##################
#WORKDIR ${TOOLS}
#RUN git clone https://gitlab.com/ezlab/busco.git && \
    #cd busco && \
    #python3 setup.py install && \
    #cd config && \
    #cp config.ini config.ini.default
    ##mkdir -p data && \
    ##cd data && \
    ##wget --output-document bacteria_odb10.tar.gz https://busco-data.ezlab.org/v4/data/lineages/bacteria_odb10.2020-03-06.tar.gz && \
    ##gunzip bacteria_odb10.tar.gz && \
    ##tar -xf bacteria_odb10.tar && \
    ##wget --output-document fungi_odb10.tar.gz https://busco-data.ezlab.org/v4/data/lineages/fungi_odb10.2020-09-10.tar.gz && \
    ##gunzip fungi_odb10.tar.gz && \
    ##tar -xf fungi_odb10.tar && \
    ##wget --output-document saccharomycetes_odb10.tar.gz https://busco-data.ezlab.org/v4/data/lineages/saccharomycetes_odb10.2020-08-05.tar.gz && \
    ##gunzip saccharomycetes_odb10.tar.gz && \
    ##tar -xf saccharomycetes_odb10.tar && \
    ##wget --output-document gammaproteobacteria_odb10.tar.gz https://busco-data.ezlab.org/v4/data/lineages/gammaproteobacteria_odb10.2020-03-06.tar.gz && \
    ##gunzip gammaproteobacteria_odb10.tar.gz && \
    ##tar -xf saccharomycetes_odb10.tar &&  \
    ##cd ../config && \
    ##cp config.ini config.ini.default

## Fix the busco config file line by line:
#WORKDIR ${TOOLS}/busco/config
#RUN sed -i 's:/ncbi-blast-2.10.1+:/tools/ncbi-blast-2.9.0+:g' config.ini && \
    #sed -i 's:/augustus/:/tools/augustus-3.3.2/:g' config.ini && \
    #sed -i 's:/usr/local/bin/:/usr/bin/:g' config.ini && \
    #sed -i 's:/home/biodocker/sepp/:/usr/local/bin:g' config.ini

##########
### KAT ##
##########
#WORKDIR ${TOOLS}
#RUN pip3 install tabulate
##RUN git clone https://github.com/TGAC/KAT.git \
##    && cd KAT \
##    && ./build_boost.sh \
##    && ./autogen.sh \
##    && ./configure \
##    && make 

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

##############
### Merqury ##
##############
WORKDIR ${TOOLS}
RUN wget https://github.com/marbl/merqury/archive/v1.3.tar.gz \
    && tar -zxvf v1.3.tar.gz
ENV MERQURY=${TOOLS}/merqury-1.3


#WORKDIR ${TOOLS}
#RUN prokka-1.14.0/bin/prokka --setupdb

### Set PATH and do some extraneous steps
ENV PATH=/:/usr/src/pteryx/pteryx:/tools/ntHits-ntHits-v0.0.1/:/tools/minimap2:/tools/racon/build/bin:/tools/ntEdit/:/tools/racon-v1.3.1/build/bin:/tools/augustus-3.3.2/bin:/tools/augustus-3.3.2/bin/scripts:/tools/ncbi-blast-2.9.0+/bin:/tools/Flye/bin:/tools/prokka-1.14.0/bin/:/tools/barrnap-0.8/bin:/tools/bbmap:/tools/bowtie2-2.3.2/:/tools/SPAdes-4.0.0-Linux/bin:/tools/Filtlong/bin:/tools/canu-1.8/Linux-amd64/bin:/tools:/usr/bin:/tools/SKESA:/tools/miniasm/:/tools/seqtk:/tools/quast-5.0.2:/tools/minigraph:/tools/SPAdes-3.13.0-Linux/bin:$PATH:/tools/mmseqs/bin/:$PATH

#RUN rm ${TOOLS}/*.gz && \
    #rm ${TOOLS}/*.zip

#################################################################################
# Workflow and homemade scripts
# Set up the working directory
ENV APP_HOME=/usr/src/pteryx
RUN mkdir -p $APP_HOME
WORKDIR $APP_HOME

# Copy necessary files
COPY setup.py setup.py
COPY pteryx pteryx
COPY setup.cfg setup.cfg

# Install build tool using pip (not uv) to ensure it's in the system Python path
RUN pip install build --break-system-packages

# Install build tool
RUN uv pip install build

# Ensure uv-installed binaries are in PATH
ENV PATH="/root/.local/bin:${PATH}"

# Build the package
RUN python3 -m build

# Install the built package
RUN uv pip install dist/*.whl

# Set Python path
ENV PYTHONPATH=${PYTHONPATH}:/usr/src/pteryx


# Set up application home
#ENV APP_HOME=/pteryx
#RUN mkdir -p $APP_HOME
#WORKDIR $APP_HOME

# Create input and output directories
RUN mkdir -p $APP_HOME/inputs $APP_HOME/outputs

