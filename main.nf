#!/usr/bin/env nextflow
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    REPO4EU/modulediscovery
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    Github : https://github.com/REPO4EU/modulediscovery
----------------------------------------------------------------------------------------
*/

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS / WORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { MODULEDISCOVERY  } from './workflows/modulediscovery'
include { PIPELINE_INITIALISATION } from './subworkflows/local/utils_nfcore_modulediscovery_pipeline'
include { PIPELINE_COMPLETION     } from './subworkflows/local/utils_nfcore_modulediscovery_pipeline'
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    NAMED WORKFLOWS FOR PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// WORKFLOW: Run main analysis pipeline depending on type of input
//
workflow REPO4EU_MODULEDISCOVERY {

    take:
    ch_seeds                // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    ch_network              // channel: [ val(meta[id,network_id]), path(network) ]
    ch_shortest_paths       // channel: [ val(meta[id,network_id]), path(shortest_paths) ]
    ch_permuted_networks    // channel: [ val(meta[id,network_id]), [path(permuted_networks)] ]

    main:

    ch_versions = Channel.empty()

    //
    // WORKFLOW: Run pipeline
    //

    MODULEDISCOVERY (
        ch_seeds,
        ch_network,
        ch_shortest_paths,
        ch_permuted_networks
    )
    ch_versions = ch_versions.mix(MODULEDISCOVERY.out.versions)

    emit:
    multiqc_report = MODULEDISCOVERY.out.multiqc_report // channel: /path/to/multiqc_report.html
    versions       = ch_versions                        // channel: [version1, version2, ...]

}
/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow {

    main:
    //
    // SUBWORKFLOW: Run initialisation tasks
    //
    PIPELINE_INITIALISATION (
        params.version,
        params.validate_params,
        params.monochrome_logs,
        args,
        params.outdir,
        params.input,
        params.seeds,
        params.network,
        params.shortest_paths,
        params.permuted_networks,
        params.id_space
    )

    //
    // WORKFLOW: Run main workflow
    //
    REPO4EU_MODULEDISCOVERY (PIPELINE_INITIALISATION.out.seeds, PIPELINE_INITIALISATION.out.network, PIPELINE_INITIALISATION.out.shortest_paths, PIPELINE_INITIALISATION.out.permuted_networks)

    //
    // SUBWORKFLOW: Run completion tasks
    //
    PIPELINE_COMPLETION (
        params.email,
        params.email_on_fail,
        params.plaintext_email,
        params.outdir,
        params.monochrome_logs,
        params.hook_url,
        REPO4EU_MODULEDISCOVERY.out.multiqc_report
    )
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
