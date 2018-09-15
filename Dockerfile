#################################################################
# Dockerfile
#
# Version:          2.0.10
# Software:         QIMERA base image
# Software Version: 20180828
# Description:      Dependencies for QIMERA
# Website:          http://www.qiagen.com/
# Tags:             NGS|RNA|Fusion|QIAGEN|Transcriptomics
# Base Image:       debian:9
# Build Cmd:        docker build -t qiagen/qimera -f Dockerfile .
#################################################################

# base image
FROM python:3

# personal info
MAINTAINER Frank Reinecke <frank.reinecke@qiagen.com>

# Set noninterative mode
ENV DEBIAN_FRONTEND noninteractive

# apt update and install global requirements
RUN apt-get clean all && \
    apt-get update && \
    apt-get upgrade -y && \
    apt-get install -y apt-utils

RUN apt-get install -y --no-install-recommends \
        automake         \
        autotools-dev    \
        build-essential  \
        bzip2            \
        ca-certificates  \
        cmake            \
        curl             \
        dpkg             \
        git              \
        grep             \
        gzip             \
        make             \
        pkg-config       \
        pmtools          \
        python           \
        sed              \
        unzip            \
        wget             \
        zip              \
        zlib1g-dev && \
        apt-get clean

# install system packages of perl modules (faster than cpanminus)
RUN apt-get install -y -q    \
     libarchive-zip-perl           \
     libdbd-sqlite3-perl           \
     libdbd-mysql-perl             \
     liblog-log4perl-perl          \
     libmoose-perl                 \
     libexcel-writer-xlsx-perl     \
     libparallel-forkmanager-perl  \
     libwww-perl                   \
     libmodule-build-perl          \
     libmethod-signatures-perl     \
     libmoosex-configfromfile-perl \
     libdatetime-perl   && \
     apt-get clean

# install miniconda
RUN echo 'export PATH=/opt/conda/bin:$PATH' > /etc/profile.d/conda.sh && \
    wget --quiet https://repo.continuum.io/miniconda/Miniconda2-4.0.5-Linux-x86_64.sh -O ~/miniconda.sh && \
    /bin/bash ~/miniconda.sh -b -p /opt/conda && \
    rm ~/miniconda.sh

RUN TINI_VERSION=`curl https://github.com/krallin/tini/releases/latest | grep -o "/v.*\"" | sed 's:^..\(.*\).$:\1:'` && \
    curl -L "https://github.com/krallin/tini/releases/download/v${TINI_VERSION}/tini_${TINI_VERSION}.deb" > tini.deb && \
    dpkg -i tini.deb && \
    rm tini.deb && \
    apt-get clean

RUN mkdir /data && mkdir /references

# Add qiauser 
RUN useradd --create-home --shell /bin/bash --user-group --uid 1000 --groups sudo qiauser && \
    echo `echo "qiauser\nqiauser\n" | passwd qiauser` && \
    chown qiauser:qiauser /data && \
    chown qiauser:qiauser /references && \
    chmod 777 -R /opt/conda/

# Add $HOME/bin to path
ENV PATH=$PATH:/opt/conda/bin
ENV PATH=$PATH:/home/qiauser/bin
ENV HOME=/home/qiauser

# Create $HOME/bin folder
RUN mkdir /home/qiauser/bin

# Add R and bioconda channel
RUN conda config --add channels r && \
    conda config --add channels bioconda && \
    conda upgrade conda

######## BEGIN INSTALLATION OF SPECIFIC TOOLS ##############
#
# install bedtools, BLAST, cutadapt, samtools, STAR, SPAdes
RUN conda install bedtools=2.25.0 cutadapt=1.9.1 samtools=1.3.1 star=2.5.2a blast=2.2.31 spades=3.10.1

# install seq-align from GitHub 
RUN git clone --recursive https://github.com/noporpoise/seq-align && \
	cd seq-align && \
	make && \
	cp bin/* /bin/

# install cpanminus
RUN curl -L http://cpanmin.us | perl - --self-upgrade

# insall deps only
RUN curl -s https://api.github.com/repos/qiaseq/QIMERA/releases/latest | \
        grep "QIMERA-" | cut -d : -f 2,3 | tr -d \" | tail -1 | \
	cpanm -q --installdeps --without-suggests --without-recommends --skip-satisfied --no-man-pages

# # install own code now
# RUN curl -s https://api.github.com/repos/qiaseq/QIMERA/releases/latest | \
#         grep "QIMERA-" | cut -d : -f 2,3 | tr -d \" | tail -1 | \
#         cpanm -q --notest && rm -rf /home/qiauser/.cpanm/work

# trimmer deps
RUN pip3 install edlib Cython

# smCounter deps
RUN wget https://bootstrap.pypa.io/get-pip.py -O /tmp/get-pip.py && python2 /tmp/get-pip.py

RUN apt-get -y install --no-install-recommends \
	python-dev \
 	zlib1g-dev \
	bzip2 \
	liblzma-dev \
	libcurl4-openssl-dev

RUN pip2 install numpy scipy pysam==0.9.0
RUN apt-get -y install r-base
RUN echo "r <- getOption('repos'); r['CRAN'] <- 'http://cran.us.r-project.org'; options(repos = r);" > ~/.Rprofile
RUN Rscript -e "install.packages('plyr')"

# QC pipeline dependencies
RUN pip3 install joblib

# clean up
RUN apt-get clean all && apt-get purge && rm -rf /var/lib/apt/lists/* /tmp/* /var/tmp/*

######### END INSTALLATION OF SPECIFIC TOOLS ##############

ENTRYPOINT []
CMD []

##################### INSTALLATION END #####################
