<h1>
  <picture>
    <source media="(prefers-color-scheme: dark)" srcset="docs/images/nf-core-diseasemodulediscovery_logo_dark.png">
    <img alt="nf-core/diseasemodulediscovery" src="docs/images/nf-core-diseasemodulediscovery_logo_light.png">
  </picture>
</h1>

[![GitHub Actions CI Status](https://github.com/nf-core/diseasemodulediscovery/actions/workflows/ci.yml/badge.svg)](https://github.com/nf-core/diseasemodulediscovery/actions/workflows/ci.yml)
[![GitHub Actions Linting Status](https://github.com/nf-core/diseasemodulediscovery/actions/workflows/linting.yml/badge.svg)](https://github.com/nf-core/diseasemodulediscovery/actions/workflows/linting.yml)[![AWS CI](https://img.shields.io/badge/CI%20tests-full%20size-FF9900?labelColor=000000&logo=Amazon%20AWS)](https://nf-co.re/diseasemodulediscovery/results)[![Cite with Zenodo](http://img.shields.io/badge/DOI-10.5281/zenodo.XXXXXXX-1073c8?labelColor=000000)](https://doi.org/10.5281/zenodo.XXXXXXX)
[![nf-test](https://img.shields.io/badge/unit_tests-nf--test-337ab7.svg)](https://www.nf-test.com)

[![Nextflow](https://img.shields.io/badge/nextflow%20DSL2-%E2%89%A524.04.2-23aa62.svg)](https://www.nextflow.io/)
[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
[![run with singularity](https://img.shields.io/badge/run%20with-singularity-1d355c.svg?labelColor=000000)](https://sylabs.io/docs/)
[![Launch on Seqera Platform](https://img.shields.io/badge/Launch%20%F0%9F%9A%80-Seqera%20Platform-%234256e7)](https://cloud.seqera.io/launch?pipeline=https://github.com/nf-core/diseasemodulediscovery)

[![Get help on Slack](http://img.shields.io/badge/slack-nf--core%20%23diseasemodulediscovery-4A154B?labelColor=000000&logo=slack)](https://nfcore.slack.com/channels/diseasemodulediscovery)[![Follow on Twitter](http://img.shields.io/badge/twitter-%40nf__core-1DA1F2?labelColor=000000&logo=twitter)](https://twitter.com/nf_core)[![Follow on Mastodon](https://img.shields.io/badge/mastodon-nf__core-6364ff?labelColor=FFFFFF&logo=mastodon)](https://mstdn.science/@nf_core)[![Watch on YouTube](http://img.shields.io/badge/youtube-nf--core-FF0000?labelColor=000000&logo=youtube)](https://www.youtube.com/c/nf-core)

## Introduction

**nf-core/diseasemodulediscovery** is a bioinformatics pipeline for network medicine hypothesis generation, designed for identifying active/disease modules. Developed and maintained by the [RePo4EU](https://repo4.eu/) consortium, it aims to characterize the molecular mechanisms of diseases by analyzing the local neighborhood of disease-associated genes or proteins (seeds) within the interactome. This approach can help identify potential drug targets for drug repurposing.

![REPO4EU/modulediscovery metro map](docs/images/REPO4EU_modulediscovery_metro_map.png)

- Module inference (all enabled by default):
  - [`DOMINO`](https://github.com/Shamir-Lab/DOMINO)
  - [`DIAMOnD`](https://github.com/dinaghiassian/DIAMOnD)
  - [`ROBUST`](https://github.com/bionetslab/robust)
  - [`ROBUST bias aware`](https://github.com/bionetslab/robust_bias_aware)
  - `first neighbors`
  - `random walk with restart`
- Visualization of the module networks ([`graph-tool`](https://graph-tool.skewed.de/), [`pyvis`](https://github.com/WestHealth/pyvis))
- Export to the network medicine web visualization tool [`Drugst.One`](https://drugst.one/)
- Annotation with biological data (targeting drugs, side effects, associated disorders, cellular localization) queried from [`NeDRexDB`](https://nedrex.net/) and conversion to [`BioPAX`](https://www.biopax.org/) format.
- Evaluation
  - Over-representation analysis ([`g:Profiler`](https://cran.r-project.org/web/packages/gprofiler2/index.html))
  - Functional coherence analysis ([`DIGEST`](https://pypi.org/project/biodigest/))
  - Network topology analysis ([`graph-tool`](https://graph-tool.skewed.de/))
  - Seed set permutation-based evaluation (enabled by `--run_seed_permutation`)
  - Network permutation-based evaluation (enabled by `--run_network_permutation`)
- Drug prioritization using the API of [`Drugst.One`](https://drugst.one/)
- Result and evaluation summary ([`MultiQC`](https://seqera.io/multiqc/))

## Usage

> [!NOTE]
> If you are new to Nextflow and nf-core, please refer to [this page](https://nf-co.re/docs/usage/installation) on how to set-up Nextflow. Make sure to [test your setup](https://nf-co.re/docs/usage/introduction#how-to-run-a-pipeline) with `-profile test` before running the workflow on actual data.

### Running the pipeline

Now, you can run the pipeline using:

```bash
nextflow run nf-core/diseasemodulediscovery \
   -profile <docker/singularity> \
   --seeds <SEED_FILE> \
   --network <NETWORK_FILE> \
   --outdir <OUTDIR>
```

This will run the pipeline based on the provided `<SEED_FILE>` and `<NETWORK_FILE>`. Results will be saved to the specified `<OUTDIR>`. Use `-profile` to set whether docker or singularity should be used for software deployment.

> [!WARNING]
> Please provide pipeline parameters via the CLI or Nextflow `-params-file` option. Custom config files including those provided by the `-c` Nextflow option can be used to provide any configuration _**except for parameters**_; see [docs](https://nf-co.re/docs/usage/getting_started/configuration#custom-configuration-files).

For more details and further functionality, please refer to the [usage documentation](https://nf-co.re/diseasemodulediscovery/usage) and the [parameter documentation](https://nf-co.re/diseasemodulediscovery/parameters).

## Pipeline output

To see the results of an example test run with a full size dataset refer to the [results](https://nf-co.re/diseasemodulediscovery/results) tab on the nf-core website pipeline page.
For more details about the output files and reports, please refer to the
[output documentation](https://nf-co.re/diseasemodulediscovery/output).

## Credits

nf-core/diseasemodulediscovery was originally written by the [RePo4EU](https://repo4.eu/) consortium.

We thank the following people for their extensive assistance in the development of this pipeline:

- [Johannes Kersting](https://github.com/JohannesKersting) (TUM)
- [Lisa Spindler](https://github.com/lspindler2509) (TUM)
- [Quirin Manz](https://github.com/quirinmanz) (TUM)
- [Quim Aguirre](https://github.com/quimaguirre) (STALICLA)
- [Chloé Bucheron](https://github.com/ChloeBubu) (University Vienna)

## Contributions and Support

If you would like to contribute to this pipeline, please see the [contributing guidelines](.github/CONTRIBUTING.md).

If you want to include an additional module identification approach, please see [this guide](docs/contributing.md).
For further information or help, don't hesitate to get in touch on the [Slack `#diseasemodulediscovery` channel](https://nfcore.slack.com/channels/diseasemodulediscovery) (you can join with [this invite](https://nf-co.re/join/slack)).

## Citations

<!-- TODO nf-core: Add citation for pipeline after first release. Uncomment lines below and update Zenodo doi and badge at the top of this file. -->
<!-- If you use nf-core/diseasemodulediscovery for your analysis, please cite it using the following doi: [10.5281/zenodo.XXXXXX](https://doi.org/10.5281/zenodo.XXXXXX) -->

<!-- TODO nf-core: Add bibliography of tools and data used in your pipeline -->

An extensive list of references for the tools used by the pipeline can be found in the [`CITATIONS.md`](CITATIONS.md) file.

You can cite the `nf-core` publication as follows:

> **The nf-core framework for community-curated bioinformatics pipelines.**
>
> Philip Ewels, Alexander Peltzer, Sven Fillinger, Harshil Patel, Johannes Alneberg, Andreas Wilm, Maxime Ulysse Garcia, Paolo Di Tommaso & Sven Nahnsen.
>
> _Nat Biotechnol._ 2020 Feb 13. doi: [10.1038/s41587-020-0439-x](https://dx.doi.org/10.1038/s41587-020-0439-x).
