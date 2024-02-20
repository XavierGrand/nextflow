# Getting Started

These instructions will get you a copy of the project up and running on your local machine for development and testing purposes.
You can follow the [building your pipeline guide](./doc/building_your_pipeline.md) to learn how to build your own pipelines.

## Prerequisites

To run nextflow on your computer you need to have `java` (>= 1.8) installed.

```sh
java --version
```

and `git`

```sh
git --version
```

To be able to run existing tools in nextflow on your computer (`src/nf_modules/` to see the list). You need to have `docker` installed.

```sh
docker run hello-world
```

Alternatively if you are on Linux, you can use `singularity`:

```sh
singularity run docker://hello-world
```

## Installing

To install nextflow on your computer simply run the following command:

```sh
git clone git@gitbio.ens-lyon.fr:LBMC/nextflow.git
cd nextflow/
src/install_nextflow.sh
```

## Running a toy RNASeq quantification pipeline

To run tests we first need to get a training set
```sh
cd data
git clone https://gitbio.ens-lyon.fr/LBMC/Hub/tiny_dataset.git
cp tiny_dataset/fastq/tiny_R1.fastq tiny_dataset/fastq/tiny2_R1.fastq
cp tiny_dataset/fastq/tiny_R2.fastq tiny_dataset/fastq/tiny2_R2.fastq
cp tiny_dataset/fastq/tiny_S.fastq tiny_dataset/fastq/tiny2_S.fastq
cd ..
```

Then to run the tests for a given tools run the following command:

```sh
./nextflow src/solution_RNASeq.nf --fastq "data/tiny_dataset/fastq/tiny2_R{1,2}.fastq.gz" --fasta "data/tiny_dataset/fasta/tiny_v2_10.fasta" --bed "data/tiny_dataset/annot/tiny.bed" -profile docker
```

## Nextflow profile

By default le `src/nextflow.config` file define 4 different profiles

- `-profile docker` each process of the pipeline will be executed within a `docker` container locally
- `-profile singularity` each process of the pipeline will be executed within a `singularity` container locally
- `-profile psmn` each process will be sent as a separate job within a `singularity` container on the PSMN
- `-profile ccin2p3` each process will be sent as a separate job within a `singularity` container on the CCIN2P3

If the containers are not found locally, they are automatically downloaded before running the process. For the PSMN and CCIN2P3, the `singularity` images are downloaded in a shared folder (`/scratch/Bio/singularity` for the PSMN, and `/sps/lbmc/common/singularity/` for the CCIN2P3)


### PSMN

When running `nextflow` on the PSMN, we recommend to use `tmux` before launching the pipeline:

```sh
tmux
./nextflow src/solution_RNASeq.nf --fastq "data/tiny_dataset/fastq/tiny2_R{1,2}.fastq.gz" --fasta "data/tiny_dataset/fasta/tiny_v2_10.fasta" --bed "data/tiny_dataset/annot/tiny.bed" -profile psmn
```

Therefore, the `nextflow` process will continue to run even if you are disconnected.
You can re-attach the `tmux` session, with the command `tmux a` (and press `ctrl` `+` `b` `+` `d` to detach the attached session).

### CCIN2P3

When runnning `nextflow` on the CCIN2P3, you cannot use `tmux`, instead you should send a *daemon* jobs which will launch the `nextflow` command.
You can edit the `src/ccin2p3.pbs` file to personalize your `nextflow` command and send it with the command:

```sh
qsub src/ccin2p3.pbs
```

