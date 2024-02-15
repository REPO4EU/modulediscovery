[![GitHub Actions CI Status](https://github.com/REPO4EU/modulediscovery/workflows/nf-core%20CI/badge.svg)](https://github.com/REPO4EU/modulediscovery/actions?query=workflow%3A%22nf-core+CI%22)
[![GitHub Actions Linting Status](https://github.com/REPO4EU/modulediscovery/workflows/nf-core%20linting/badge.svg)](https://github.com/REPO4EU/modulediscovery/actions?query=workflow%3A%22nf-core+linting%22)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A523.04.0-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Nextflow Tower](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Nextflow%20Tower-%234256e7)](https://tower.nf/launch?pipeline=https://github.com/REPO4EU/modulediscovery)

## Introduction

**REPO4EU/modulediscovery** is a bioinformatics pipeline that ...

<!-- TODO nf-core:
   Complete this sentence with a 2-3 sentence summary of what types of data the pipeline ingests, a brief overview of the
   major pipeline sections and the types of output it produces. You're giving an overview to someone new
   to nf-core here, in 15-20 seconds. For an example, see https://github.com/nf-core/rnaseq/blob/master/README.md#introduction
-->

<!-- TODO nf-core: Include a figure that guides the user through the major workflow steps. Many nf-core
     workflows use the "tube map" design for that. See https://nf-co.re/docs/contributing/design_guidelines#examples for examples.   -->
<!-- TODO nf-core: Fill in short bullet-pointed list of the default steps in the pipeline -->

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

Clone the repository:

```bash
git clone https://github.com/REPO4EU/modulediscovery.git
cd modulediscovery
```

Install required software e.g. via conda (docker not included):

```bash
conda env create -n modulediscovery -f environment.yaml
conda activate modulediscovery
```

The pipeline should be run from outside of the code repository since nextflow, by default, will write into the execution directory.

Run with test data (DIAMOnD example data):

```bash
nextflow run <PATH_TO_REPO>/modulediscovery/main.nf \
   -profile <docker/singularity>,test \
   --outdir <OUTDIR>
```

Now, you can run the pipeline using:

<!-- TODO nf-core: update the following command to include all required parameters for a minimal example -->

```bash
nextflow run <PATH_TO_REPO>/modulediscovery/main.nf \
   -profile <docker/singularity> \
   --input <seed_file> \
   --network <network_file> \
   --outdir <OUTDIR>
```

Show all parameter options:

```bash
nextflow run <PATH_TO_REPO>/modulediscovery/main.nf --help
```

If you want to contribute to the pipeline, it is useful to set up pre-commit for code linting and quality checks:

```bash
pre-commit install
```

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_;
> see [docs](https://nf-co.re/usage/configuration#custom-configuration-files).

## Including a new active module detection tool

1. Create a new branch for your tool.
2. Add a function to the bin/graph_tool_parser.py script for preparing the tool-specific network input format. The script is built around the [graph-tool](https://graph-tool.skewed.de/) Python package. An example is the safe_diamond() function, which saves a simple edge list in CSV format. Add the function as an option in the save() function and a command line option for `--format` in parse_args(). The output file name has to include the option specified with `--format` since nextflow uses this pattern to check whether the output file was successfully generated. The script expects a .gt file as input. Run the pipeline with the "test" profile to generate a .gt example file in `<OUTDIR>/graphtoolparser`, which you can use to test the parsing function by executing the parsing script directly via the command line.
3. Create a module for the tool. (Example with comments: `modules/local/diamond/main.nf` and `modules/local/domino/`)
4. Create a subworkflow wrapping the tool together with the input parser. (Example with comments: `subworkflows/local/gt_diamond/main.nf` and `subworkflows/local/gt_domino/main.nf`)
5. Include the subworkflow in the `workflows/modulediscovery.nf` file. Again, DIAMOnD and DOMINO are included as examples.
6. Test checks locally:
   1. Run tests via, e.g., `nextflow run main.nf -profile singularity,test --outdir results`.
   2. Run `nf-core lint`.
   3. Check your code style. This will automatically happen before you commit, if you use pre-commit, which can be set up with: `pre-commit install`. After each commit, it will automatically check your code style and fix it where possible. If changes were made, you have to commit again.
7. Create a pull request against the dev branch.

### Further information

- [FAQ sheet](https://docs.google.com/document/d/1WgBIFrrcxFKN0I-zJbuS7PUCmyCLPTWx6xAHg1zi4FA/edit?usp=sharing)
- [Workflow schema](https://docs.google.com/drawings/d/1X7U79dAZaeRdGdIsXoEKw74MNqjxCHq3RuNASBYCiB4/edit?usp=sharing)

## Credits

REPO4EU/modulediscovery was originally written by REPO4EU.

We thank the following people for their extensive assistance in the development of this pipeline:

<!-- TODO nf-core: If applicable, make list of people who have also contributed -->

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use REPO4EU/modulediscovery for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

This pipeline uses code and infrastructure developed and maintained by the [nf-core](https://nf-co.re) community, reused here under the [MIT license](https://github.com/nf-core/tools/blob/master/LICENSE).

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
