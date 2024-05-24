/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Loaded from modules/local/
//
include { GRAPHTOOLPARSER         } from '../modules/local/graphtoolparser/main'
include { GT2TSV as GT2TSV_Modules} from '../modules/local/gt2tsv/main'
include { GT2TSV as GT2TSV_Network} from '../modules/local/gt2tsv/main'
include { ADDHEADER               } from '../modules/local/addheader/main'
include { DIGEST                  } from '../modules/local/digest/main'


//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { GT_BIOPAX         } from '../subworkflows/local/gt_biopax/main'
include { NETWORKEXPANSION  } from '../subworkflows/local/networkexpansion/main'
include { PERMUTATION       } from '../subworkflows/local/permutation/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//                                                //Evaluation
include { GPROFILER2_GOST        } from '../modules/nf-core/gprofiler2/gost/main'
include { MULTIQC                } from '../modules/nf-core/multiqc/main'
include { paramsSummaryMap       } from 'plugin/nf-validation'
include { paramsSummaryMultiqc   } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { softwareVersionsToYAML } from '../subworkflows/nf-core/utils_nfcore_pipeline'
include { methodsDescriptionText } from '../subworkflows/local/utils_nfcore_modulediscovery_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow MODULEDISCOVERY {


    take:
    ch_seeds // channel: samplesheet read in from --input
    ch_network // channel: network file read in from --network

    main:

    // Params
    id_space = Channel.value(params.id_space)
    validate_online = Channel.value(params.validate_online)

    // Channels
    ch_versions = Channel.empty()
    ch_multiqc_files = Channel.empty()


    // Brach channel, so, GRAPHTOOLPARSER runs only for supported network formats, which are not already .gt files
    ch_network_type = ch_network.branch {
        gt: it.extension == "gt"
        parse: true
    }

    // Run network parser for non .gt networks, supported by graph-tool
    GRAPHTOOLPARSER(ch_network_type.parse, 'gt')
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)

    // Mix into one .gt format channel
    ch_network_gt = GRAPHTOOLPARSER.out.network.collect().mix(ch_network_type.gt).first()


    // Network expansion tools
    NETWORKEXPANSION(ch_seeds, ch_network_gt)
    ch_modules = NETWORKEXPANSION.out.modules
    ch_versions = ch_versions.mix(NETWORKEXPANSION.out.versions)

    // Permutation based evaluation
    PERMUTATION(ch_seeds, ch_modules, ch_network_gt)
    ch_versions = ch_versions.mix(PERMUTATION.out.versions)


    // Annotation and BIOPAX conversion
    if(!params.skip_annotation){
        GT_BIOPAX(ch_modules, id_space, validate_online)
        ch_versions = ch_versions.mix(GT_BIOPAX.out.versions)
    }

    GT2TSV_Modules(ch_modules)
    GT2TSV_Network(ch_network_gt.flatten().map{ it -> [ [ id: it.baseName ], it ] })
    ADDHEADER(ch_seeds, "gene_id")
    ch_nodes = GT2TSV_Modules.out
    ch_nodes = ch_nodes.mix(ADDHEADER.out)
    // Evaluation
    if(!params.skip_gprofiler){

        GPROFILER2_GOST (
            ch_nodes,
            [],
            GT2TSV_Network.out.map{it[1]}.first()
        )
        ch_versions = ch_versions.mix(GPROFILER2_GOST.out.versions)
    }

    if(!params.skip_digest){
        DIGEST (ch_nodes, id_space, ch_network_gt, id_space)
        ch_versions = ch_versions.mix(DIGEST.out.versions)
    }

    // Collate and save software versions
    softwareVersionsToYAML(ch_versions)
        .collectFile(
            storeDir: "${params.outdir}/pipeline_info",
            name: 'nf_core_pipeline_software_mqc_versions.yml',
            sort: true,
            newLine: true
        ).set { ch_collated_versions }

    //
    // MODULE: MultiQC
    //
    ch_multiqc_config        = Channel.fromPath(
        "$projectDir/assets/multiqc_config.yml", checkIfExists: true)
    ch_multiqc_custom_config = params.multiqc_config ?
        Channel.fromPath(params.multiqc_config, checkIfExists: true) :
        Channel.empty()
    ch_multiqc_logo          = params.multiqc_logo ?
        Channel.fromPath(params.multiqc_logo, checkIfExists: true) :
        Channel.empty()

    summary_params      = paramsSummaryMap(
        workflow, parameters_schema: "nextflow_schema.json")
    ch_workflow_summary = Channel.value(paramsSummaryMultiqc(summary_params))

    ch_multiqc_custom_methods_description = params.multiqc_methods_description ?
        file(params.multiqc_methods_description, checkIfExists: true) :
        file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)
    ch_methods_description                = Channel.value(
        methodsDescriptionText(ch_multiqc_custom_methods_description))

    ch_multiqc_files = ch_multiqc_files.mix(
        ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_collated_versions)
    ch_multiqc_files = ch_multiqc_files.mix(
        ch_methods_description.collectFile(
            name: 'methods_description_mqc.yaml',
            sort: true
        )
    )

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )

    emit:
    multiqc_report = MULTIQC.out.report.toList() // channel: /path/to/multiqc_report.html
    versions       = ch_versions                 // channel: [ path(versions.yml) ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
