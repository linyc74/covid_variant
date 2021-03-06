FROM continuumio/miniconda3:4.9.2

RUN apt-get install -y \
    libtbb2 \
 && apt-get clean

RUN conda install --channel bioconda \
    bowtie2=2.4.1 \
    samtools=1.11 \
    bcftools=1.11 \
 && conda clean --all --yes

RUN pip install --no-cache-dir \
    pandas==1.2.4 \
    biopython==1.79 \
    ngslite==1.2.1 \
    cutadapt==3.3

COPY ./covid_variant/* /covid_variant/covid_variant/
COPY ./__main__.py /covid_variant/
WORKDIR /
ENTRYPOINT ["python", "covid_variant"]
