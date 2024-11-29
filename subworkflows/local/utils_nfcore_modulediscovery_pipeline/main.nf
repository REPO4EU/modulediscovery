//
// Subworkflow with functionality specific to the REPO4EU/modulediscovery pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFVALIDATION_PLUGIN } from '../../nf-core/utils_nfvalidation_plugin'
include { paramsSummaryMap          } from 'plugin/nf-validation'
include { fromSamplesheet           } from 'plugin/nf-validation'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { dashedLine                } from '../../nf-core/utils_nfcore_pipeline'
include { nfCoreLogo                } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { workflowCitation          } from '../../nf-core/utils_nfcore_pipeline'

/*
========================================================================================
    SUBWORKFLOW TO INITIALISE PIPELINE
========================================================================================
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    help              // boolean: Display help text
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to sample sheet
    seeds             //  string: Path(s) to seed file(s)
    network           //  string: Path(s) to network file(s)
    shortest_paths    //  string: Path to shortest paths file

    main:

    ch_versions = Channel.empty()

    //
    // Print version and exit if required and dump pipeline parameters to JSON file
    //
    UTILS_NEXTFLOW_PIPELINE (
        version,
        true,
        outdir,
        workflow.profile.tokenize(',').intersect(['conda', 'mamba']).size() >= 1
    )

    //
    // Validate parameters and generate parameter summary to stdout
    //
    pre_help_text = nfCoreLogo(monochrome_logs)
    post_help_text = '\n' + workflowCitation() + '\n' + dashedLine(monochrome_logs)
    def String workflow_command = "nextflow run ${workflow.manifest.name} -profile <docker/singularity/.../institute> --seeds <seed_file> --network <network_file> --outdir <OUTDIR>"
    UTILS_NFVALIDATION_PLUGIN (
        help,
        workflow_command,
        pre_help_text,
        post_help_text,
        validate_params,
        "nextflow_schema.json"
    )

    //
    // Check config provided to the pipeline
    //
    UTILS_NFCORE_PIPELINE (
        nextflow_cli_args
    )

    ch_seeds = Channel.empty()          // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    ch_network = Channel.empty()        // channel: [ val(meta[id,network_id]), path(network) ]
    ch_shortest_paths = Channel.empty() // channel: [ val(meta[id,network_id]), path(shortest_paths) ]

    seed_param_set = (params.seeds != null)
    network_param_set = (params.network != null)
    shortest_paths_param_set = (params.shortest_paths != null)

    if(params.input){

        // Check sample sheet

        // channel: [ path(seeds), path(network) ]
        ch_input = Channel
            .fromSamplesheet("input")
            .map{seeds, network, shortest_paths ->
                if((seeds.size()==0) ^ seed_param_set ){
                    error("Seed genes have to specified through either the sample sheet OR the --seeds parameter")
                }
                if((network.size()==0) ^ network_param_set){
                    error("Networks have to specified through either the sample sheet OR the --network parameter")
                }
                if(!(shortest_paths.size()==0) && shortest_paths_param_set ){
                    error("Shortest paths have to specified through either the sample sheet OR the --shortest_path parameter")
                }
                if(!(network.size()==0) && shortest_paths_param_set ){
                    error("If the network is set via the sample sheet, shortest_paths must also be set via the sample sheet")
                }
                if(!(shortest_paths.size()==0) && network_param_set ){
                    error("If the shortest_paths is set via the sample sheet, the network must also be set via the sample sheet")
                }
                [seeds, network, shortest_paths]
            }

        if (seed_param_set && network_param_set) {

            error("You need to specify either a sample sheet (--input) OR the seeds (--seeds) and network (--network) files")

        } else if (!seed_param_set && !network_param_set) {

            log.info("Creating network and seeds channels based on tuples in the sample sheet")

            ch_network = ch_input
                .map{ it -> [it[1], it[2]]}
                .map{ network, sp ->
                    network: [ [ id: network.baseName, network_id: network.baseName ], network, sp ]
                }
                .unique()

            ch_seeds = ch_input
                .map{ it ->
                    seeds = it[0]
                    network = it[1]
                    network_id = network.baseName
                    [ [ id: seeds.baseName + "." + network_id, seeds_id: seeds.baseName + "." + network_id, network_id: network_id ] , seeds ]
                }

        } else if (seed_param_set && !network_param_set) {

            log.info("Creating network channel based on the sample sheet and seeds channel based on the seeds parameter")

            ch_network = ch_input
                .map{ it -> it[1]}
                .map{ [ [ id: it.baseName, network_id: it.baseName ], it ] }

            ch_seeds = Channel
                .fromPath(params.seeds.split(',').flatten(), checkIfExists: true)
                .combine(ch_network.map{meta, network -> meta.network_id})
                .map{seeds, network_id ->
                    [ [ id: seeds.baseName + "." + network_id, seeds_id: seeds.baseName + "." + network_id, network_id: network_id ] , seeds ]
                }

        } else if (!seed_param_set && network_param_set) {

            log.info("Creating network channel based on the network parameter and seeds channel based on the sample sheet")

            ch_network = Channel
                .fromPath(params.network.split(',').flatten(), checkIfExists: true)
                .map{ it -> [ [ id: it.baseName, network_id: it.baseName ], it ] }

            ch_seeds = ch_input
                .map{ it -> it[0]}
                .combine(ch_network.map{meta, network -> meta.network_id})
                .map{seeds, network_id ->
                    [ [ id: seeds.baseName + "." + network_id, seeds_id: seeds.baseName + "." + network_id, network_id: network_id ] , seeds ]
                }

        }


    } else if (seed_param_set && network_param_set){

        log.info("Creating network and seeds channels based on the combination of all seed and network files provided")

        ch_network = Channel
            .fromPath(params.network.split(',').flatten(), checkIfExists: true)
            .map{ it -> [ [ id: it.baseName, network_id: it.baseName ], it ] }

        ch_seeds = Channel
            .fromPath(params.seeds.split(',').flatten(), checkIfExists: true)
            .combine(ch_network.map{meta, network -> meta.network_id})
            .map{seeds, network_id ->
                [ [ id: seeds.baseName + "." + network_id, seeds_id: seeds.baseName + "." + network_id, network_id: network_id ] , seeds ]
            }

    } else {
        error("You need to specify either a sample sheet (--input) or the seeds (--seeds) and network (--network) files")
    }

    // check if IDs are unique
    ch_network.map{ meta, network, sp -> [meta.id] }
        .collect()
        .subscribe { list ->
            def unique = list.size() == list.toSet().size()
            if (!unique) { error("IDs in ch_network are not unique.") }
        }
    ch_seeds.map{ meta, seeds -> [meta.id] }
        .collect()
        .subscribe { list ->
            def unique = list.size() == list.toSet().size()
            if (!unique) { error("IDs in ch_seeds are not unique.") }
        }

    // separate shortest paths
    ch_shortest_paths = ch_network.map{meta, network, sp -> [meta, sp.size() > 0 ? sp : file("${projectDir}/assets/NO_FILE", checkIfExists: true)]}
    ch_network = ch_network.map{meta, network, sp -> [meta, network]}

    emit:
    versions    = ch_versions
    seeds       = ch_seeds             // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    network     = ch_network           // channel: [ val(meta[id,network_id]), path(network) ]
    shortest_paths = ch_shortest_paths // channel: [ val(meta[id,network_id]), path(shortest_paths) ]
}

/*
========================================================================================
    SUBWORKFLOW FOR PIPELINE COMPLETION
========================================================================================
*/

workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report

    main:

    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")

    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(summary_params, email, email_on_fail, plaintext_email, outdir, monochrome_logs, multiqc_report.toList())
        }

        completionSummary(monochrome_logs)

        if (hook_url) {
            imNotification(summary_params, hook_url)
        }
    }

    workflow.onError {
        log.error "Pipeline failed. Please refer to troubleshooting docs: https://nf-co.re/docs/usage/troubleshooting"
    }
}

/*
========================================================================================
    FUNCTIONS
========================================================================================
*/

//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ it.single_end }.unique().size == 1
    if (!endedness_ok) {
        error("Please check input samplesheet -> Multiple runs of a sample must be of the same datatype i.e. single-end or paired-end: ${metas[0].id}")
    }

    return [ metas[0], fastqs ]
}

//
// Generate methods description for MultiQC
//
def toolCitationText() {
    // TODO nf-core: Optionally add in-text citation tools to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "Tool (Foo et al. 2023)" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def citation_text = [
            "Tools used in the workflow included:",
            "FastQC (Andrews 2010),",
            "MultiQC (Ewels et al. 2016)",
            "."
        ].join(' ').trim()

    return citation_text
}

def toolBibliographyText() {
    // TODO nf-core: Optionally add bibliographic entries to this list.
    // Can use ternary operators to dynamically construct based conditions, e.g. params["run_xyz"] ? "<li>Author (2023) Pub name, Journal, DOI</li>" : "",
    // Uncomment function in methodsDescriptionText to render in MultiQC report
    def reference_text = [
            "<li>Andrews S, (2010) FastQC, URL: https://www.bioinformatics.babraham.ac.uk/projects/fastqc/).</li>",
            "<li>Ewels, P., Magnusson, M., Lundin, S., & Käller, M. (2016). MultiQC: summarize analysis results for multiple tools and samples in a single report. Bioinformatics , 32(19), 3047–3048. doi: /10.1093/bioinformatics/btw354</li>"
        ].join(' ').trim()

    return reference_text
}

def methodsDescriptionText(mqc_methods_yaml) {
    // Convert  to a named map so can be used as with familar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        String[] manifest_doi = meta.manifest_map.doi.tokenize(",")
        for (String doi_ref: manifest_doi) temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        meta["doi_text"] = temp_doi_ref.substring(0, temp_doi_ref.length() - 2)
    } else meta["doi_text"] = ""
    meta["nodoi_text"] = meta.manifest_map.doi ? "" : "<li>If available, make sure to update the text to include the Zenodo DOI of version of the pipeline used. </li>"

    // Tool references
    meta["tool_citations"] = ""
    meta["tool_bibliography"] = ""

    // TODO nf-core: Only uncomment below if logic in toolCitationText/toolBibliographyText has been filled!
    // meta["tool_citations"] = toolCitationText().replaceAll(", \\.", ".").replaceAll("\\. \\.", ".").replaceAll(", \\.", ".")
    // meta["tool_bibliography"] = toolBibliographyText()


    def methods_text = mqc_methods_yaml.text

    def engine =  new groovy.text.SimpleTemplateEngine()
    def description_html = engine.createTemplate(methods_text).make(meta)

    return description_html.toString()
}
