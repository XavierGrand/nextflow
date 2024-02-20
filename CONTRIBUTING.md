# Contributing

When contributing to this repository, please first discuss the change you wish to make via issues,
email, or on the [ENS-Bioinfo channel](https://matrix.to/#/#ens-bioinfo:matrix.org) before making a change. 

## Forking

In git, the [action of forking](https://git-scm.com/book/en/v2/GitHub-Contributing-to-a-Project) means that you are going to make your own private copy of a repository. You can then write modifications in your project, and if they are of interest for the source repository create a merge request (here [LBMC/nextflow](https://gitbio.ens-lyon.fr/LBMC/nextflow)). Merge requests are sent to the source repository to ask the maintainers to integrate modifications.

![merge request button](./doc/img/merge_request.png)

## Project organization

The `LBMC/nextflow` project is structured as follows:
- all the code is in the `src/` folder
- scripts downloading external tools should download them in the `bin/` folder
- all the documentation (including this file) can be found int he `doc/` folder
- the `data` and `results` folders contain the data and results of your pipelines and are ignored by `git`

## Code structure

The `src/` folder is where we want to save the pipeline (`.nf`) scripts. This folder also contains
- the `src/install_nextflow.sh` to install the nextflow executable at the root of the project.
- some pipelines examples (like the one build during the nf_pratical)
- the `src/nextflow.config` global configuration file which contains the `docker`, `singularity`, `psmn` and `ccin2p3` profiles.
- the `src/nf_modules` folder contains per tools `main.nf` modules with predefined process that users can import in their projects with the [DSL2](https://www.nextflow.io/docs/latest/dsl2.html)

But also some hidden folders that users don't need to see when building their pipeline:
- the `src/.docker_modules` contains the recipes for the `docker` containers used in the `src/nf_modules/<tool_names>/main.nf` files
- the `src/.singularity_in2p3` and `src/.singularity_psmn` are symbolic links to the shared folder where the singularity images are downloaded on the PSMN and CCIN2P3 

# Proposing a new tool

Each tool named `<tool_name>` must have two dedicated folders:

- [`src/nf_modules/<tool_name>`](./src/nf_modules/fastp/) where users can find `.nf` files to include
- [`src/.docker_modules/<tool_name>/<version_number>`](./src/.docker_modules/fastp/0.20.1/) where we have the [`Dockerfile`](./src/.docker_modules/fastp/0.20.1/Dockerfile) to construct the container used in the `main.nf` file

## `src/nf_module` guide lines

We are going to take the [`fastp`, `nf_module`](./src/nf_modules/fastp/) as an example.

The [`src/nf_modules/<tool_name>`](./src/nf_modules/fastp/) should contain a [`main.nf`](./src/nf_modules/fastp/main.nf) file that describe at least one process using `<tool_name>`

### container informations

The first two lines of [`main.nf`](./src/nf_modules/fastp/main.nf) should define two variables
```Groovy
version = "0.20.1"
container_url = "lbmc/fastp:${version}"
```

we can then use the `container_url` definition in each `process` in the `container` attribute.
In addition to the `container` directive, each `process` should have one of the following `label` attributes (defined in the `src/nextflow.config` file)
- `big_mem_mono_cpus`
- `big_mem_multi_cpus`
- `small_mem_mono_cpus`
- `small_mem_multi_cpus`

```Groovy
process fastp {
  container = "${container_url}"
  label = "big_mem_multi_cpus"
  ...
}
```

### process options

Before each process, you should declare at least two `params.` variables:
- A `params.<process_name>` defaulting to `""` (empty string) to allow user to add more command line option to your process without rewriting the process definition
- A `params.<process_name>_out` defaulting to `""` (empty string) that define the `results/` subfolder where the process output should be copied if the user wants to save the process output

```Groovy
params.fastp = ""
params.fastp_out = ""
process fastp {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  if (params.fastp_out != "") {
    publishDir "results/${params.fastp_out}", mode: 'copy'
  }
  ...
  script:
"""
fastp --thread ${task.cpus} \
${params.fastp} \
...
"""
}
```

The user can then change the value of these variables:
- from the command line `--fastp "--trim_head1=10"``
- with the `include` command within their pipeline: `include { fastq } from "nf_modules/fastq/main" addParams(fastq_out: "QC/fastq/")
- by defining the variable within their pipeline: `params.fastq_out = "QC/fastq/"

### `input` and `output` format

You should always use `tuple` for input and output channel format with at least:
- a `val` containing variable(s) related to the item
- a `path` for the file(s) that you want to process

for example:

```Groovy
process fastp {
  container = "${container_url}"
  label "big_mem_multi_cpus"
  tag "$file_id"
  if (params.fastp_out != "") {
    publishDir "results/${params.fastp_out}", mode: 'copy'
  }

  input:
  tuple val(file_id), path(reads)

  output:
    tuple val(file_id), path("*.fastq.gz"), emit: fastq
    tuple val(file_id), path("*.html"), emit: html
    tuple val(file_id), path("*.json"), emit: report
...
```

Here `file_id` can be anything from a simple identifier to a list of several variables.
In which case the first item of the List should be usable as a file prefix.
So you have to keep that in mind if you want to use it to define output file names (you can test for that with `file_id instanceof List`).
In some case, the `file_id` may be a Map to have a cleaner access to the `file_id` content by explicit keywords.

If you want to use information within the `file_id` to name outputs in your `script` section, you can use the following snipet:

```Groovy
  script:
    switch(file_id) {
    case {it instanceof List}:
      file_prefix = file_id[0]
    break
    case {it instanceof Map}:
      file_prefix = file_id.values()[0]
    break
    default:
      file_prefix = file_id
    break
  }
```

and use the `file_prefix` variable.

This also means that channel emitting `path` item should be transformed with at least the following map function:

```Groovy
.map { it -> [it.simpleName, it]}
```

for example

```Groovy
channel
  .fromPath( params.fasta )
  .ifEmpty { error "Cannot find any fasta files matching: ${params.fasta}" }
  .map { it -> [it.simpleName, it]}
  .set { fasta_files }
```


The rationale behind taking a `file_id` and emitting the same `file_id` is to facilitate complex channel operations in pipelines without having to rewrite the `process` blocks.

### dealing with paired-end and single-end data

When oppening fastq files with `channel.fromFilePairs( params.fastq )`, item in the channel have the following shape:

```Groovy
[file_id, [read_1_file, read_2_file]]
```

To make this call more generic, we can use the `size: -1` option, and accept arbitrary number of associated fastq files:

```Groovy
channel.fromFilePairs( params.fastq, size: -1 )
```

will thus give `[file_id, [read_1_file, read_2_file]]` for paired-end data and `[file_id, [read_1_file]]` for single-end data


You can the use tests on `read.size()` to define conditional `script` block:

```Groovy
...
  script:
  if (file_id instanceof List){
    file_prefix = file_id[0]
  } else {
    file_prefix = file_id
  }
  if (reads.size() == 2)
  """
  fastp --thread ${task.cpus} \
    ${params.fastp} \
    --in1 ${reads[0]} \
    --in2 ${reads[1]} \
    --out1 ${file_prefix}_R1_trim.fastq.gz \
    --out2 ${file_prefix}_R2_trim.fastq.gz \
    --html ${file_prefix}.html \
    --json ${file_prefix}_fastp.json \
    --report_title ${file_prefix}
  """
  else
  """
  fastp --thread ${task.cpus} \
    ${params.fastp} \
    --in1 ${reads[0]} \
    --out1 ${file_prefix}_trim.fastq.gz \
    --html ${file_prefix}.html \
    --json ${file_prefix}_fastp.json \
    --report_title ${file_prefix}
  """
...
```

### Complex processes

Sometime you want to do write complex processes, for example for `fastp` we want to have predefine `fastp` process for different protocols, order of adapter trimming and reads clipping.
We can then use the fact that `process` or named `workflow` can be interchangeably imported with th [DSL2](https://www.nextflow.io/docs/latest/dsl2.html#workflow-composition).

With the following example, the user can simply include the `fastp` step without knowing that it's a named `workflow` instead of a `process`.
By specifying the `params.fastp_protocol`, the `fastp` step will transparently switch betwen the different `fastp` `process`es.
Here `fastp_default` or `fastp_accel_1splus`, and other protocols can be added later, pipeline will be able to handle these new protocols by simply updating from the `upstream` repository without changing their codes.

```Groovy
params.fastp_protocol = ""
workflow fastp {
  take:
    fastq

  main:
  switch(params.fastp_protocol) {
    case "accel_1splus":
      fastp_accel_1splus(fastq)
      fastp_accel_1splus.out.fastq.set{res_fastq}
      fastp_accel_1splus.out.report.set{res_report}
    break;
    default:
      fastp_default(fastq)
      fastp_default.out.fastq.set{res_fastq}
      fastp_default.out.report.set{res_report}
    break;
  }
  emit:
    fastq = res_fastq
    report = res_report
}
```

## `src/.docker_modules` guide lines

We are going to take the [`fastp`, `.docker_modules`](./src/.docker_module/fastp/0.20.1/) as an example.

The [`src/.docker_modules/<tool_name>/<version_number>`](./src/nf_modules/fastp/0.20.1/) should contain a [`Dockerfile`](./src/.docker_module/fastp/0.20.1/Dockerfile) and a [`docker_init.sh`](./src/.docker_module/fastp/0.20.1/docker_init.sh).

### `Dockerfile`

The [`Dockerfile`](./src/.docker_module/fastp/0.20.1/Dockerfile) shoud contains a `docker` recipe to build a image with `<tool_name>` installed in a system-wide binary folder (`/bin`, `/usr/local/bin/`, etc).
Therefore, your scripts are easily accessible from within the container.

This recipe should have:

- an easily changeable `<version_number>` to be able to update the corresponding image to a newer version of the tool
- the `ps` executable (package `procps` in debian)
- a default `bash` command (`CMD ["bash"]`)

### `docker_init.sh`

The [`docker_init.sh`](./src/.docker_module/fastp/0.20.1/docker_init.sh) script is a small sh script with the following content:

```sh
#!/bin/sh
docker pull lbmc/fastp:0.20.1
docker build src/.docker_modules/fastp/0.20.1 -t 'lbmc/fastp:0.20.1'
docker push lbmc/fastp:0.20.1
```

We want to be able to execute the `src/.docker_module/fastp/0.20.1/docker_init.sh` from the root of the project to :

- try to download the corresponding container if it exists on the [Docker Hub](https://hub.docker.com/repository/docker/lbmc/)
- if not build the container from the correspondig [`Dockerfile`](./src/.docker_module/fastp/0.20.1/Dockerfile) and with the same name as the name we would get from the `docker pull` command
- push the container on the [Docker Hub](https://hub.docker.com/repository/docker/lbmc/) (only [laurent.modolo@ens-lyon.fr](mailto:laurent.modolo@ens-lyon.fr) can do this step for the group **lbmc**)

