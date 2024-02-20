# nextflow pipeline

This repository is a template and a library repository to help you build nextflow pipeline.
You can fork this repository to build your own pipeline.

## Getting the last updates

To get the last commits from this repository into your fork use the following commands:

For the first time:
```sh
git remote add upstream git@gitbio.ens-lyon.fr:LBMC/nextflow.git
git pull upstream master
```

Then to make an update:
```sh
git pull upstream master
git merge upstream/master
```

## Getting Started

These instructions will get you a copy of the project as a template when you want to build your own pipeline.

[you can follow them here.](doc/getting_started.md)

## Building your pipeline

You can follow the [building your pipeline guide](./doc/building_your_pipeline.md) for a gentle introduction to `nextflow` and taking advantage of this template to build your pipelines.

## Existing Nextflow pipeline

Before starting a new project, you can check if someone else didn’t already to the work !
- [on the nextflow project page](./doc/nf_projects.md)
- [on the nf-core project](https://nf-co.re/pipelines)

## Contributing

If you want to add more tools to this project, please read the [CONTRIBUTING.md](CONTRIBUTING.md).

## Authors

* **Laurent Modolo** - *Initial work*

See also the list of [contributors](https://gitbio.ens-lyon.fr/pipelines/nextflow/graphs/master) who participated in this project.

## License

This project is licensed under the CeCiLL License- see the [LICENSE](LICENSE) file for details
