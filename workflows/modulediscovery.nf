/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PRINT PARAMS SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { paramsSummaryLog; paramsSummaryMap } from 'plugin/nf-validation'

def logo = NfcoreTemplate.logo(workflow, params.monochrome_logs)
def citation = '\n' + WorkflowMain.citation(workflow) + '\n'
def summary_params = paramsSummaryMap(workflow)

// Print parameter summary log to screen
log.info logo + paramsSummaryLog(workflow) + citation

WorkflowModulediscovery.initialise(params, log)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    PARAMS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
ch_seeds = Channel.fromPath(params.input, checkIfExists: true)
ch_network = Channel.fromPath(params.network, checkIfExists: true).first()

diamond_n = Channel.value(params.diamond_n)
diamond_alpha = Channel.value(params.diamond_alpha)

rwr_scaling = Channel.value(params.rwr_scaling).map{it ? 1 : 0}
rwr_symmetrical = Channel.value(params.rwr_symmetrical).map{it ? 1 : 0}
rwr_r = Channel.value(params.rwr_r)

id_space = Channel.value(params.id_space)
validate_online = Channel.value(params.validate_online)


/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    CONFIG FILES
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

ch_multiqc_config          = Channel.fromPath("$projectDir/assets/multiqc_config.yml", checkIfExists: true)
ch_multiqc_custom_config   = params.multiqc_config ? Channel.fromPath( params.multiqc_config, checkIfExists: true ) : Channel.empty()
ch_multiqc_logo            = params.multiqc_logo   ? Channel.fromPath( params.multiqc_logo, checkIfExists: true ) : Channel.empty()
ch_multiqc_custom_methods_description = params.multiqc_methods_description ? file(params.multiqc_methods_description, checkIfExists: true) : file("$projectDir/assets/methods_description_template.yml", checkIfExists: true)

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT LOCAL MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// SUBWORKFLOW: Consisting of a mix of local and nf-core/modules
//
include { GRAPHTOOLPARSER   } from '../modules/local/graphtoolparser/main'
include { INPUT_CHECK       } from '../subworkflows/local/input_check'
include { GT2TSV as GT2TSV_Modules} from '../modules/local/gt2tsv/main'    //Evaluation
include {GT2TSV as GT2TSV_Network} from '../modules/local/gt2tsv/main'    //Evaluation
include {ADDHEADER} from '../modules/local/addheader/main' //Evaluation

include { GT_DIAMOND        } from '../subworkflows/local/gt_diamond'
include { GT_DOMINO         } from '../subworkflows/local/gt_domino'
include { GT_ROBUST         } from '../subworkflows/local/gt_robust'
include { GT_FIRSTNEIGHBOR  } from '../subworkflows/local/gt_firstneighbor'
include { GT_RWR            } from '../subworkflows/local/gt_rwr'

include { GT_BIOPAX         } from '../subworkflows/local/gt_biopax/main'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT NF-CORE MODULES/SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

//
// MODULE: Installed directly from nf-core/modules
//
include { MULTIQC                     } from '../modules/nf-core/multiqc/main'
include { CUSTOM_DUMPSOFTWAREVERSIONS } from '../modules/nf-core/custom/dumpsoftwareversions/main'                                                   //Evaluation
include { GPROFILER2_GOST } from '../modules/nf-core/gprofiler2/gost/main'  //Evaluation

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    RUN MAIN WORKFLOW
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

// Info required for completion email and summary
def multiqc_report = []

workflow MODULEDISCOVERY {

    ch_versions = Channel.empty()
    ch_modules = Channel.empty()
    ch_nodes = Channel.empty() //Evaluation


    // Brach channel, so, GRAPHTOOLPARSER runs only for supported network formats, which are not already .gt files
    ch_network_type = ch_network.branch {
        gt: it.extension == "gt"
        parse: true
    }

    // Run network parser for non .gt networks, supported by graph-tool
    GRAPHTOOLPARSER(ch_network_type.parse, 'gt')
    ch_versions = ch_versions.mix(GRAPHTOOLPARSER.out.versions)

    // Mix into one .gt format channel
    ch_network_gt = GRAPHTOOLPARSER.out.network.collect().mix(ch_network_type.gt).collect()


    // Network expansion tools
    GT_DIAMOND(ch_seeds, ch_network_gt, diamond_n, diamond_alpha)
    ch_versions = ch_versions.mix(GT_DIAMOND.out.versions)
    ch_modules = ch_modules.mix(GT_DIAMOND.out.module)

    GT_DOMINO(ch_seeds, ch_network_gt)
    ch_versions = ch_versions.mix(GT_DOMINO.out.versions)
    ch_modules = ch_modules.mix(GT_DOMINO.out.module)

    GT_ROBUST(ch_seeds, ch_network_gt)
    ch_versions = ch_versions.mix(GT_ROBUST.out.versions)
    ch_modules = ch_modules.mix(GT_ROBUST.out.module)

    GT_FIRSTNEIGHBOR(ch_seeds, ch_network_gt)
    ch_versions = ch_versions.mix(GT_FIRSTNEIGHBOR.out.versions)
    ch_modules = ch_modules.mix(GT_FIRSTNEIGHBOR.out.module)

    GT_RWR(ch_seeds, ch_network_gt, rwr_scaling, rwr_symmetrical, rwr_r)
    ch_versions = ch_versions.mix(GT_RWR.out.versions)
    ch_modules = ch_modules.mix(GT_RWR.out.module)


    // Annotation and BIOPAX conversion
    if(!params.skip_annotation){
        GT_BIOPAX(ch_modules, id_space, validate_online)
        ch_versions = ch_versions.mix(GT_BIOPAX.out.versions)
    }


    // Evaluation
    GT2TSV_Modules(ch_modules)
    GT2TSV_Network(ch_network_gt)
    ADDHEADER(ch_seeds, "gene_id")

    ch_nodes = ch_nodes.mix(GT2TSV_Modules.out)
    ch_nodes = ch_nodes.mix(ADDHEADER.out)

    ch_gprofiler_input = ch_nodes.map{[[id: it.baseName],it]}
    GPROFILER2_GOST (
        ch_gprofiler_input,
        [],
        GT2TSV_Network.out
    )
    ch_versions = ch_versions.mix(GPROFILER2_GOST.out.versions)


    // Collect software versions
    CUSTOM_DUMPSOFTWAREVERSIONS (
        ch_versions.unique().collectFile(name: 'collated_versions.yml')
    )

    //
    // MODULE: MultiQC
    //
    workflow_summary    = WorkflowModulediscovery.paramsSummaryMultiqc(workflow, summary_params)
    ch_workflow_summary = Channel.value(workflow_summary)

    methods_description    = WorkflowModulediscovery.methodsDescriptionText(workflow, ch_multiqc_custom_methods_description, params)
    ch_methods_description = Channel.value(methods_description)

    ch_multiqc_files = Channel.empty()
    ch_multiqc_files = ch_multiqc_files.mix(ch_workflow_summary.collectFile(name: 'workflow_summary_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(ch_methods_description.collectFile(name: 'methods_description_mqc.yaml'))
    ch_multiqc_files = ch_multiqc_files.mix(CUSTOM_DUMPSOFTWAREVERSIONS.out.mqc_yml.collect())

    MULTIQC (
        ch_multiqc_files.collect(),
        ch_multiqc_config.toList(),
        ch_multiqc_custom_config.toList(),
        ch_multiqc_logo.toList()
    )
    multiqc_report = MULTIQC.out.report.toList()
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    COMPLETION EMAIL AND SUMMARY
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow.onComplete {
    if (params.email || params.email_on_fail) {
        NfcoreTemplate.email(workflow, params, summary_params, projectDir, log, multiqc_report)
    }
    NfcoreTemplate.dump_parameters(workflow, params)
    NfcoreTemplate.summary(workflow, params, log)
    if (params.hook_url) {
        NfcoreTemplate.IM_notification(workflow, params, summary_params, projectDir, log)
    }
}

workflow.onError {
    if (workflow.errorReport.contains("Process requirement exceeds available memory")) {
        println("🛑 Default resources exceed availability 🛑 ")
        println("💡 See here on how to configure pipeline: https://nf-co.re/docs/usage/configuration#tuning-workflow-resources 💡")
    }
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    THE END
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/
