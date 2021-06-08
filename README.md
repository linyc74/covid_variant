# Covid Variant

**Detect COVID-19 variant using Illumina paired-end sequencing**

### Usage

    git clone https://github.com/linyc74/covid_variant.git

    python covid_variant -1 read1.fq.gz -2 read2.fq.gz

### Reference Sequence and Variants

- WT COVID-19 genome: [`NC_045512.2.gb`](https://www.ncbi.nlm.nih.gov/nuccore/1798174254)
- [SARS-CoV-2 Variant Classifications and Definitions](https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html)

### Output

An example shown in `stdout`

    Spike protein mutations: 69del, 70del, 143del, 144del, Y145D, N501Y, A570D, D614G, P681H, T716I, A942S, S982A, D1118H
    Match: B.1.1.7 [United Kingdom]

### Dependency

Download and install Anaconda on either Mac or Linux. For Windows users, Windows Subsystem for Linux works as well.

Once Anaconda is set up, install the following packages in the terminal:

    pip install ngslite

    conda install -c bioconda \
        trim_galore \
        bowtie2 \
        samtools \
        bcftools
