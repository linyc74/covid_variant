# Covid Variant

**Detect COVID-19 variant using Illumina sequencing**

### Usage

    git clone https://github.com/linyc74/covid_variant.git

Paired-end

    python covid_variant -1 read1.fq.gz -2 read2.fq.gz

Unpaired

    python covid_variant -1 read.fq.gz

### Reference Sequence and Variants

- WT COVID-19 genome: [`NC_045512.2.gb`](https://www.ncbi.nlm.nih.gov/nuccore/1798174254)
- [SARS-CoV-2 Variant Classifications and Definitions](https://www.cdc.gov/coronavirus/2019-ncov/variants/variant-info.html)

|Name     |Spike Protein Substitutions                                              |First Detected            |
|---------|-------------------------------------------------------------------------|--------------------------|
|B.1.1.7  |69del, 70del, 144del, N501Y, A570D, D614G, P681H, T716I, S982A, D1118H   |United Kingdom            |
|B.1.351  |D80A, D215G, 241del, 242del, 243del, K417N, E484K, N501Y, D614G, A701V   |South Africa              |
|B.1.427  |L452R, D614G                                                             |United States-(California)|
|B.1.429  |S13I, W152C, L452R, D614G                                                |United States-(California)|
|B.1.617  |L452R, E484Q, D614G                                                      |India – February 2021     |
|B.1.617.1|G142D, E154K, L452R, E484Q, D614G, P681R, Q1071H                         |India – December 2020     |
|B.1.617.2|T19R, 156del, 157del, R158G, L452R, T478K, D614G, P681R, D950N           |India – December 2020     |
|B.1.617.3|T19R, G142D, L452R, E484Q, D614G, P681R, D950N                           |India – October 2020      |
|P.1      |L18F, T20N, P26S, D138Y, R190S, K417T, E484K, N501Y, D614G, H655Y, T1027I|Japan/Brazil              |
|P.2      |E484K, D614G, V1176F                                                     |Brazil – April 2020       |

### Output

An example shown in `stdout`

    Spike protein mutations: 69del, 70del, 143del, 144del, Y145D, N501Y, A570D, D614G, P681H, T716I, A942S, S982A, D1118H
    Match: B.1.1.7 [United Kingdom]

### Dependency

Download and install Anaconda on either Mac or Linux. Windows Subsystem for Linux (WSL) works as well.

Once Anaconda is set up, install the following packages in the terminal:

    pip install pandas biopython ngslite cutadapt
    conda install -c bioconda fastqc bowtie2 samtools bcftools trim-galore

For some users, installing `trim-galore` might cause some trouble.
Please refer to [TrimGalore](https://github.com/FelixKrueger/TrimGalore) for more details.

### Docker

The following section is the detailed instruction for running docker container.

Pull the latest docker image from Docker Hub

    docker pull linyc74/covid-variant

Check if the container works by printing the help message

    docker run linyc74/covid-variant --help

To run docker, the folder containing fastq files must be mounted on the docker container.
For example, if the fastq files are located in the `data` folder:

    data/read1.fq.gz
    data/read2.fq.gz

Then, the `data` folder must be mounted as `--volume $(pwd)/data:/data`, where `$(pwd)` provides the absolute path.
The docker command would be:

    docker run \
    --volume $(pwd)/data:/data \
    linyc74/covid-variant \ 
    -1 /data/read1.fq.gz \
    -2 /data/read2.fq.gz

If the output directory is also needed, simply add new volumes

    docker run \
    --volume $(pwd)/data:/data \
    --volume $(pwd)/outdir:/outdir \
    linyc74/covid-variant \ 
    -1 /data/read1.fq.gz \
    -2 /data/read2.fq.gz \
    -o /outdir
