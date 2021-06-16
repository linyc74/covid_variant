## Docker Container

The following is the detailed instruction for running docker container.

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
