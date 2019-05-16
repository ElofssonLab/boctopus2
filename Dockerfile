FROM nanjiang/common-ubuntu
LABEL maintainer "Nanjiang Shu" nanjiang.shu@nbis.se
LABEL version "1.0"

# Install Miniconda2
RUN apt-get update && apt-get -y dist-upgrade && \
    apt-get install -y --no-install-recommends bzip2 curl ca-certificates language-pack-en  fontconfig vim && \
    curl https://repo.anaconda.com/miniconda/Miniconda2-latest-Linux-x86_64.sh -O && \
    bash  Miniconda2-latest-Linux-x86_64.sh -bf -p /opt/miniconda2/ && \
    rm Miniconda2-latest-Linux-x86_64.sh

# Add Conda to PATH
ENV PATH="/opt/miniconda2/bin:${PATH}"

ENV LC_ALL en_US.UTF-8
ENV LC_LANG en_US.UTF-8
SHELL ["/bin/bash", "-c"]

WORKDIR /app

RUN mkdir -p /data

# Set up the Conda environment
COPY conda_env.yml .
RUN conda config --add channels conda-forge && \
    conda env update -n base -f conda_env.yml && \
    conda clean --all

# Install HHblits packages and Legacy blast 2.2.26
RUN curl -LJO https://github.com/ElofssonLab/boctopus2/raw/master/src/app/hhsuite-2.0.16-linux-x86_64.tar.gz  && \
        tar -xzf hhsuite-2.0.16-linux-x86_64.tar.gz  && \
        rm -f  hhsuite-2.0.16-linux-x86_64.tar.gz && \
    curl ftp://ftp.ncbi.nih.gov/blast/executables/legacy.NOTSUPPORTED/2.2.26/blast-2.2.26-x64-linux.tar.gz -O && \
        tar -xzf blast-2.2.26-x64-linux.tar.gz  && \
        rm -f  blast-2.2.26-x64-linux.tar.gz 

# Add the workflow files
ADD src ./boctopus2

# Set config file for boctopus2
RUN printf "HHBLITS_PATH: /app/hhsuite-2.0.16-linux-x86_64\nHHBLITS_DB_PATH: /data/hhsuite/uniprot20_2013_03/uniprot20_2013_03\nblastpath: /app/blast-2.2.26/bin\nrpath: /opt/miniconda2/bin/\n" > boctopus2/config.yml


ENV USER_DIRS "/app"

CMD ["/bin/bash" ]
