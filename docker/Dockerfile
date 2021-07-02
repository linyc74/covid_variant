FROM continuumio/miniconda3:4.9.2

RUN apt-get install -y \
    libtbb2 \
 && conda install --channel bioconda \
    fastqc=0.11.9 \
    bowtie2=2.4.1 \
    samtools=1.9 \
    bcftools=1.9 \
 && conda clean --all --yes \
 && pip install --no-cache-dir \
    pandas==1.2.4 \
    biopython==1.79 \
    ngslite==1.2.1 \
    cutadapt==3.3

COPY ./covid_variant /covid_variant
WORKDIR /
ENTRYPOINT ["python", "covid_variant"]
