# Building your own pipeline

The goal of this guide is to walk you through the Nextflow pipeline building process you will learn:

1. How to use this [git repository (LBMC/nextflow)](https://gitbio.ens-lyon.fr/LBMC/nextflow) as a template for your project.
2. The basis of [Nextflow](https://www.nextflow.io/) the pipeline manager that we use at the lab.
3. How to build a simple pipeline for the transcript-level quantification of RNASeq data
4. How to run the exact same pipeline on a computing center ([PSMN](http://www.ens-lyon.fr/PSMN/doku.php))

This guide assumes that you followed the [Git basis, training course](https://gitbio.ens-lyon.fr/LBMC/hub/formations/git_basis).

# Initialize your own project

You are going to build a pipeline for you or your team. So the first step is to create your own project.

## Forking

Instead of reinventing the wheel, you can use the [LBMC/nextflow](https://gitbio.ens-lyon.fr/LBMC/nextflow) as a template.
To easily do so, go to the [LBMC/nextflow](https://gitbio.ens-lyon.fr/LBMC/nextflow) repository and click on the [**fork**](https://gitbio.ens-lyon.fr/LBMC/nextflow/forks/new) button (you need to log-in).

![fork button](./img/fork.png)

In git, the [action of forking](https://git-scm.com/book/en/v2/GitHub-Contributing-to-a-Project) means that you are going to make your own private copy of a repository.
This repository will keep a link with the original [LBMC/nextflow](https://gitbio.ens-lyon.fr/LBMC/nextflow) project from which you will be able to

- [get updates](https://gitbio.ens-lyon.fr/LBMC/nextflow#getting-the-last-updates) `LBMC/nextflow` from the repository
- propose update (see [contributing guide](https://gitbio.ens-lyon.fr/LBMC/nextflow/-/blob/master/CONTRIBUTING.md#forking))


## Project organization

This project (and yours) follows the [guide of good practices for the LBMC](http://www.ens-lyon.fr/LBMC/intranet/services-communs/pole-bioinformatique/ressources/good_practice_LBMC)

You are now on the main page of your fork of the [LBMC/nextflow](https://gitbio.ens-lyon.fr/LBMC/nextflow). You can explore this project, all the codes in it is under the CeCILL licence (in the [LICENCE](https://gitbio.ens-lyon.fr/LBMC/nextflow/blob/master/LICENSE) file).

The [README.md](https://gitbio.ens-lyon.fr/LBMC/nextflow/blob/master/README.md) file contains instructions to run your pipeline and test its installation.

The [CONTRIBUTING.md](https://gitbio.ens-lyon.fr/LBMC/nextflow/blob/master/CONTRIBUTING.md) file contains guidelines if you want to contribute to the [LBMC/nextflow](https://gitbio.ens-lyon.fr/LBMC/nextflow).

The [data](https://gitbio.ens-lyon.fr/LBMC/nextflow/tree/master/data) folder will be the place where you store the raw data for your analysis.
The [results](https://gitbio.ens-lyon.fr/LBMC/nextflow/tree/master/results) folder will be the place where you store the results of your analysis.

**The content of `data` and `results` folders should never be saved on git.**

The [doc](https://gitbio.ens-lyon.fr/LBMC/nextflow/tree/master/doc) folder contains the documentation and this guide.

And most interestingly for you, the [src](https://gitbio.ens-lyon.fr/LBMC/nextflow/tree/master/src) contains code to wrap tools. This folder contains one visible subdirectories `nf_modules` some pipeline examples and other hidden folders and files.

# Nextflow pipeline

A pipeline is a succession of [**process**](https://www.nextflow.io/docs/latest/process.html#process-page). Each `process` has data input(s) and optional data output(s). Data flows are modeled as [**channels**](https://www.nextflow.io/docs/latest/channel.html).

## Processes

Here is an example of **process**:

```Groovy
process sample_fasta {
  input:
  path fasta

  output:
  path "sample.fasta", emit: fasta_sample

  script:
"""
head ${fasta} > sample.fasta
"""
}
```

We have the process `sample_fasta` that takes fasta `path` input and as output a fasta `path`. The `process` task itself is defined in the `script:` block and within `"""`.

```Groovy
input:
path fasta
```

When we zoom on the `input:` block, we see that we define a variable `fasta` of type `path`.
This means that the `sample_fasta` `process` is going to get a flux of fasta file(s).
Nextflow is going to write a file named as the content of the variable `fasta` in the root of the folder where `script:` is executed.

```Groovy
output:
path "sample.fasta", emit: fasta_sample
```

At the end of the script, a file named `sample.fasta` is found in the root the folder where `script:` is executed and will be emitted as `fasta_sample`.

Using the WebIDE of Gitlab, create a file `src/fasta_sampler.nf`
![webide](./img/webide.png)

The first line that you need to add is

```Groovy
nextflow.enable.dsl=2
```

Then add the `sample_fastq` process and commit it to your repository.


## Workflow

In Nexflow, `process` blocks are chained together within a `workflow` block.
For the time being, we only have one `process` so `workflow` may look like an unnecessary complication, but keep in mind that we want to be able to write complex bioinformatic pipeline.

```Groovy
workflow {
  sample_fasta(fasta_file)
}
```

Like `process` blocks `workflow` can take some inputs: `fasta_files`
and transmit this input to `process`es

```Groovy
  sample_fasta(fasta_file)
```

Add the definition of the `workflow` to the `src/fasta_sampler.nf` file and commit it to your repository.

## Channels

Why bother with `channel`s? In the above example, the advantages of `channel`s are not really clear. We could have just given the `fasta` file to the `workflow`. But what if we have many fasta files to process? What if we have sub processes to run on each of the sampled fasta files? Nextflow can easily deal with these problems with the help of `channel`s.

> **Channels** are streams of items that are emitted by a source and consumed by a process. A process with a `channel` as input will be run on every item send through the `channel`.

```Groovy
channel
  .fromPath( "data/tiny_dataset/fasta/*.fasta" )
  .set { fasta_file }
```

Here we defined the `channel`, `fasta_file`, that is going to send every fasta file from the folder `data/tiny_dataset/fasta/` into the process that takes it as input.

Add the definition of the `channel`, above the `workflow` block, to the `src/fasta_sampler.nf` file and commit it to your repository.

## Run your pipeline locally

After writing this first pipeline, you may want to test it. To do that, first clone your repository.
After following the [Git basis, training course](https://gitbio.ens-lyon.fr/LBMC/hub/formations/git_basis), you should have an up-to-date `ssh` configuration to connect to the `gitbio.ens-lyon.fr` git server.

You can run the following commands to download your project on your computer:

```sh
git clone git@gitbio.ens-lyon.fr:<usr_name>/nextflow.git
cd nextflow
src/install_nextflow.sh
```

We also need data to test our pipeline:

```sh
cd data
git clone git@gitbio.ens-lyon.fr:LBMC/hub/tiny_dataset.git
cd ..
```

We can run our pipeline with the following command:

```sh
./nextflow src/fasta_sampler.nf
```


## Getting your results

Our pipeline seems to work but we don’t know where is the `sample.fasta`. To get results out of a `process`, we need to tell nextflow to write it somewhere (we may don’t need to get every intermediate file in our results).

To do that we need to add the following line before the `input:` section:

```Groovy
publishDir "results/sampling/", mode: 'copy'
```

Every file described in the `output:` section will be copied from nextflow to the folder `results/sampling/`.

Add this to your `src/fasta_sampler.nf` file with the WebIDE and commit to your repository.
Pull your modifications locally with the command:

```sh
git pull origin master
```

You can run your pipeline again and check the content of the folder `results/sampling`.

## Fasta everywhere

We ran our pipeline on one fasta file. How would nextflow handle 100 of them? To test that we need to duplicate the `tiny_v2.fasta` file: 

```sh
for i in {1..100}
do
cp data/tiny_dataset/fasta/tiny_v2.fasta data/tiny_dataset/fasta/tiny_v2_${i}.fasta
done
```

You can run your pipeline again and check the content of the folder `results/sampling`.

Every `fasta_sampler` process write a `sample.fasta` file. We need to make the name of the output file dependent of the name of the input file.

```Groovy
output:
path "*_sample.fasta", emit: fasta_sample

  script:
"""
head ${fasta} > ${fasta.simpleName}_sample.fasta
"""
```

Add this to your `src/fasta_sampler.nf` file with the WebIDE and commit it to your repository before pulling your modifications locally.
You can run your pipeline again and check the content of the folder `results/sampling`.

Congratulations you built your first, one step, nextflow pipeline !


# Build your own RNASeq pipeline

In this section you are going to build your own pipeline for RNASeq analysis from the code available in the `src/nf_modules` folder.

Open the WebIDE and create a `src/RNASeq.nf` file.

The first line that we are going to add is

```Groovy
nextflow.enable.dsl=2
```

## fastp 

The first step of the pipeline is to remove any Illumina adaptors left in your read files and to trim your reads by quality.

The [LBMC/nextflow](https://gitbio.ens-lyon.fr/LBMC/nextflow) template provide you with many tools, for which you can find a predefined `process` block.
You can find a list of these tools in the [`src/nf_modules`](./src/nf_modules) folder.
You can also ask for a new tool by creating a [new issue for it](https://gitbio.ens-lyon.fr/LBMC/nextflow/-/issues/new) in the [LBMC/nextflow](https://gitbio.ens-lyon.fr/LBMC/nextflow) project.

We are going to include the [`src/nf_modules/fastp/main.nf`](./src/nf_modules/fastp/main.nf) in our `src/RNASeq.nf` file

```Groovy
include { fastp } from "./nf_modules/fastp/main.nf"
```
The `./nf_modules/fastp/main.nf` is relative to the `src/RNASeq.nf` file, this is why we don’t include the `src/` part of the path.

With this line we can call the `fastp` block in our future `workflow` without having to write it !
If we check the content of the file [`src/nf_modules/fastp/main.nf`](./src/nf_modules/fastp/main.nf), we can see that by including `fastp`, we are including a sub-`workflow` (we will come back on this object latter). Sub-`workflow` can be used like `process`es.

This `sub-workflow` takes a `fastq` `channel`. We need to make one:

```Groovy
channel
  .fromFilePairs( "data/tiny_dataset/fastq/*_R{1,2}.fastq", size: -1)
  .set { fastq_files }
```

The `.fromFilePairs()` function creates a `channel` of pairs of fastq files. Therefore, the items emitted by the `fastq_files` channel are going to be pairs of fastq for paired-end data.

The option `size: -1` allows for arbitrary numbers of associated files. Therefore, we can use the same `channel` creation for single-end data.

We can now include the `workflow` definition, passing the `fastq_files` `channel` to `fastp` to our `src/RNASeq.nf` file.

```Groovy
workflow {
  fastp(fastq_files)
}
```

You can commit your `src/RNASeq.nf` file, `pull` your modification locally and run your pipeline with the command:

```Groovy
./nextflow src/RNASeq.nf
```

What is happening ?

## Nextflow `-profile`

Nextflow tells you the following error: `fastp: command not found`. You haven’t `fastp` installed on your computer.

Tools installation can be a tedious process and reinstalling old version of those tools to reproduce old analyses can be very difficult.
Containers technologies like [Docker](https://www.docker.com/) or [Singularity](https://sylabs.io/singularity/) create small virtual environments where we can install software in a given version with all it’s dependencies. This environment can be saved, and shared, to have access to this exact working version of the software.

> Why two different systems ?

> Docker is easy to use and can be installed on Windows / MacOS / GNU/Linux but need admin rights.
> Singularity can only be used on GNU/Linux but don’t need admin rights, and can be used on shared environment.

The [LBMC/nextflow](https://gitbio.ens-lyon.fr/LBMC/nextflow) template provides you with [4 different `-profile`s to run your pipeline](https://gitbio.ens-lyon.fr/LBMC/nextflow/-/blob/master/doc/getting_started.md#nextflow-profile).

Profiles are defined in the [`src/nextflow.config`](./src/nextflow.config), which is the default configuration file for your pipeline (you don’t have to edit this file).

To run the pipeline locally you can use the profile `singularity` or `docker`

```Groovy
./nextflow src/RNASeq.nf -profile singularity
```

The `fastp`, `singularity` or `docker`, image is downloaded automatically and the fastq files are processed.

## Pipeline `--` arguments

We have defined the fastq files path within our `src/RNASeq.nf` file.
But what if we want to share our pipeline with someone who doesn’t want to analyze the `tiny_dataset` and but other fastq.
We can define a variable instead of fixing the path.

```Groovy
params.fastq = "data/fastq/*_{1,2}.fastq"
channel
  .fromFilePairs( params.fastq, size: -1)
  .set { fastq_files }
```

We declare a variable that contains the path of the fastq file to look for. The advantage of using `params.fastq` is that the option `--fastq` is now a parameter of your pipeline.

Thus, you can call your pipeline with the `--fastq` option.

You can commit your `src/RNASeq.nf` file, `pull` your modification locally and run your pipeline with the command:

```sh
./nextflow src/RNASeq.nf -profile singularity --fastq "data/tiny_dataset/fastq/*_R{1,2}.fastq"
```

We can also add the following line:

```Groovy
log.info "fastq files: ${params.fastq}"
```

This line simply displays the value of the variable

## BEDtools

We need the sequences of the transcripts that need to be quantified. We are going to extract these sequences from the reference `data/tiny_dataset/fasta/tiny_v2.fasta` with the `bed` file annotation `data/tiny_dataset/annot/tiny.bed`.

You can include the `fasta_from_bed` `process` from the [src/nf_modules/bedtools/main.nf](https://gitbio.ens-lyon.fr/LBMC/nextflow/blob/master/src/nf_modules/bedtools/main.nf) file to your `src/RNASeq.nf` file.

You need to be able to input a `fasta_files` `channel` and a `bed_files` `channel`.

```Groovy
log.info "fasta file : ${params.fasta}"
log.info "bed file : ${params.bed}"

channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .map { it -> [it.simpleName, it]}
  .set { fasta_files }
channel
  .fromPath( params.bed )
  .ifEmpty { error "Cannot find any bed files matching: ${params.bed}" }
  .map { it -> [it.simpleName, it]}
  .set { bed_files }
```

We introduce 2 new directives:
- `.ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }` to throw an error if the path of the file is not right
- `.map { it -> [it.simpleName, it]}` to transform our `channel` to a format compatible with the [`CONTRIBUTING`](../CONTRIBUTING.md) rules. Item, in the `channel` have the following shape [file_id, [file]], like the ones emited by the `.fromFilePairs(..., size: -1)` function.

We can add the `fastq_from_bed` step to our `workflow`

```Groovy
workflow {
  fastp(fastq_files)
  fasta_from_bed(fasta_files, bed_files)
}
```

Commit your work and test your pipeline with the following command:

```sh
./nextflow src/RNASeq.nf -profile singularity --fastq "data/tiny_dataset/fastq/*_R{1,2}.fastq" --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" --bed "data/tiny_dataset/annot/tiny.bed"
```

## Kallisto

Kallisto run in two steps: the indexation of the reference and the quantification of the transcript on this index.

You can include two `process`es with the following syntax:

```Groovy
include { index_fasta; mapping_fastq } from './nf_modules/kallisto/main.nf'
```

The `index_fasta` process needs to take as input the output of your `fasta_from_bed` `process`, which has the shape `[fasta_id, [fasta_file]]`.

The input of your `mapping_fastq` `process` needs to take as input and the output of your `index_fasta` `process` and the `fastp` `process`, of shape `[index_id, [index_file]]`, and `[fastq_id, [fastq_r1_file, fastq_r2_file]]`.

The output of a `process` is accessible through `<process_name>.out`.
In the cases where we have an `emit: <channel_name>` we can access the corresponding channel with`<process_name>.out.<channel_name>`

```Groovy
workflow {
  fastp(fastq_files)
  fasta_from_bed(fasta_files, bed_files)
  index_fasta(fasta_from_bed.out.fasta)
  mapping_fastq(index_fasta.out.index.collect(), fastp.out.fastq)
}
```

Commit your work and test your pipeline.

## Returning results

By default none of the `process` defined in `src/nf_modules` use the `publishDir` instruction.
You can specify their `publishDir` directory by specifying the :

```Groovy
params.<process_name>_out = "path"
```

Where "path" will describe a path within the `results` folder

Therefore you can either:

- call your pipeline with the following parameter `--mapping_fastq_out "quantification/"`
- add the following line to your `src/RNASeq.nf` file to get the output of the `mapping_fastq` process:

```Groovy
include { index_fasta; mapping_fastq } from './nf_modules/kallisto/main.nf' addParams(mapping_fastq_out: "quantification/")
```

Commit your work and test your pipeline.
You now have a RNASeq analysis pipeline that can run locally with Docker or Singularity!

## Bonus

A file `report.html` is created for each run with the detail of your pipeline execution.
You can use the `-resume` option to be able to save into cache the process results (the in a `work/` folder).

# Run your RNASeq pipeline on the PSMN

First you need to connect to the PSMN:

```sh
login@allo-psmn
```
Then once connected to `allo-psmn`, you can connect to `e5-2667v4comp1`:

```sh
login@m6142comp2
```

## Set your environment

Create and go to your `scratch` folder:

```sh
mkdir -p /scratch/Bio/<login>
cd /scratch/Bio/<login>
```

Then you need to clone your pipeline and get the data:

```sh
git clone https://gitbio.ens-lyon.fr/<usr_name>/nextflow.git
cd nextflow/data
git clone https://gitbio.ens-lyon.fr/LBMC/hub/tiny_dataset.git
cd ..
```

## Run nextflow

As we don’t want nextflow to be killed in case of disconnection, we start by launching `tmux`. In case of disconnection, you can restore your session with the command `tmux a` and close one with `ctr + b + d`

```sh
tmux
src/install_nextflow.sh
./nextflow src/RNASeq.nf -profile psmn --fastq "data/tiny_dataset/fastq/*_R{1,2}.fastq" --fasta "data/tiny_dataset/fasta/tiny_v2.fasta" --bed "data/tiny_dataset/annot/tiny.bed"
```

You just ran your pipeline on the PSMN!
