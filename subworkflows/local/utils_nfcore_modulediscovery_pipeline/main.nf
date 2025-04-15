//
// Subworkflow with functionality specific to the REPO4EU/modulediscovery pipeline
//

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    IMPORT FUNCTIONS / MODULES / SUBWORKFLOWS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

include { UTILS_NFSCHEMA_PLUGIN     } from '../../nf-core/utils_nfschema_plugin'
include { paramsSummaryMap          } from 'plugin/nf-schema'
include { samplesheetToList         } from 'plugin/nf-schema'
include { completionEmail           } from '../../nf-core/utils_nfcore_pipeline'
include { completionSummary         } from '../../nf-core/utils_nfcore_pipeline'
include { imNotification            } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NFCORE_PIPELINE     } from '../../nf-core/utils_nfcore_pipeline'
include { UTILS_NEXTFLOW_PIPELINE   } from '../../nf-core/utils_nextflow_pipeline'
include { logColours                } from '../../nf-core/utils_nfcore_pipeline'

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW TO INITIALISE PIPELINE
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

workflow PIPELINE_INITIALISATION {

    take:
    version           // boolean: Display version and exit
    validate_params   // boolean: Boolean whether to validate parameters against the schema at runtime
    monochrome_logs   // boolean: Do not use coloured log outputs
    nextflow_cli_args //   array: List of positional nextflow CLI args
    outdir            //  string: The output directory where the results will be saved
    input             //  string: Path to sample sheet
    seeds             //  string: Path(s) to seed file(s)
    network           //  string: Path(s) to network file(s)
    shortest_paths    //  string: Path to shortest paths file
    permuted_networks //  string: Path to folder(s) with permuted network files
    id_space          //  string: ID space to use for prepared networks

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
    UTILS_NFSCHEMA_PLUGIN (
        workflow,
        validate_params,
        null
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
    ch_permuted_networks = Channel.empty() // channel: [ val(meta[id,network_id]), [path(permuted_network)] ]

    seed_param_set = (params.seeds != null)
    network_param_set = (params.network != null)
    shortest_paths_param_set = (params.shortest_paths != null)
    permuted_networks_param_set = (params.permuted_networks != null)

    // prepare network channel, if parameter is set
    if(network_param_set){
        ch_network = Channel.fromList(params.network.split(',').flatten())
            .map{network -> mapPreparedNetwork(network, params.id_space)}
            .map{ it -> [ [ id: it.baseName, network_id: it.baseName ], it ] }
    }

    if(params.input){

        //
        // Create channel from input file provided through params.input
        //

        // channel: [ path(seeds), path(network), path(shortest_paths), path(permuted_networks) ]
        ch_input = Channel
            .fromList(samplesheetToList(params.input, "${projectDir}/assets/schema_input.json"))
            .map{seeds, network, shortest_paths, permuted_networks ->
                if((seeds.size()==0) ^ seed_param_set ){
                    error("Seed genes have to specified through either the sample sheet OR the --seeds parameter")
                }
                if((network.size()==0) ^ network_param_set){
                    error("Networks have to specified through either the sample sheet OR the --network parameter")
                }
                if(!(shortest_paths.size()==0) && shortest_paths_param_set ){
                    error("Shortest paths have to specified through either the sample sheet OR the --shortest_path parameter")
                }
                if(!(permuted_networks.size()==0) && permuted_networks_param_set ){
                    error("Precomputed network permutations have to specified through either the sample sheet OR the --permuted_networks parameter")
                }
                if(!(network.size()==0) && (shortest_paths_param_set || permuted_networks_param_set) ){
                    error("If the network is set via the sample sheet, shortest_paths or permuted_networks must also be set via the sample sheet")
                }
                if((! shortest_paths.size()==0 || ! permuted_networks.size()==0) && network_param_set ){
                    error("If the shortest_paths or permuted_networks are set via the sample sheet, the network must also be set via the sample sheet")
                }
                [seeds, network, shortest_paths, permuted_networks]
            }

        // prepare network channel, if parameter is not set
        if (!network_param_set){
            ch_network = ch_input
                .map{ it -> [it[1], it[2], it[3]]}
                .map{ network, sp, permuted_networks ->
                    [ mapPreparedNetwork(network, params.id_space), sp, permuted_networks ]
                }
                .map{ network, sp, permuted_networks ->
                    [ [ id: network.baseName, network_id: network.baseName ], network, sp, permuted_networks ]
                }
                .unique()
        }

        if (seed_param_set && network_param_set) {

            error("You need to specify either a sample sheet (--input) OR the seeds (--seeds) and network (--network) files")

        } else if (!seed_param_set && !network_param_set) {

            log.info("Creating network and seeds channels based on tuples in the sample sheet")

            ch_seeds = ch_input
                .map{ it ->
                    seeds = it[0]
                    network = it[1]
                    network_id = mapPreparedNetwork(network, params.id_space).baseName
                    [ [ id: seeds.baseName + "." + network_id, seeds_id: seeds.baseName, network_id: network_id ] , seeds ]
                }

        } else if (seed_param_set && !network_param_set) {

            log.info("Creating network channel based on the sample sheet and seeds channel based on the seeds parameter")

            ch_seeds = Channel
                .fromPath(params.seeds.split(',').flatten(), checkIfExists: true)
                .combine(ch_network.map{meta, network, sp, permuted_networks -> meta.network_id})
                .map{seeds, network_id ->
                    [ [ id: seeds.baseName + "." + network_id, seeds_id: seeds.baseName, network_id: network_id ] , seeds ]
                }

        } else if (!seed_param_set && network_param_set) {

            log.info("Creating network channel based on the network parameter and seeds channel based on the sample sheet")

            ch_seeds = ch_input
                .map{ it -> it[0]}
                .combine(ch_network.map{meta, network -> meta.network_id})
                .map{seeds, network_id ->
                    [ [ id: seeds.baseName + "." + network_id, seeds_id: seeds.baseName, network_id: network_id ] , seeds ]
                }

            // Add sp files, if provided (currently does not check if the number of the shortest paths matches the number of the networks and does not work with missing values)
            if(shortest_paths_param_set){
                ch_network = ch_network.merge(
                    Channel
                    .fromPath(params.shortest_paths.split(',').flatten())
                )
            } else{
                ch_network = ch_network.map{meta, network -> [meta, network, file("${projectDir}/assets/NO_FILE", checkIfExists: true)]}
            }

            // Add permuted network folders, if provided (currently does not check if the number of the shortest paths matches the number of the networks and does not work with missing values)
            if(permuted_networks_param_set){
                ch_network = ch_network.merge(
                    Channel
                    .fromPath(params.permuted_networks.split(',').flatten())
                )
            } else{
                ch_network = ch_network.map{meta, network, sp -> [meta, network, sp, []]}
            }

        }


    } else if (seed_param_set && network_param_set){

        log.info("Creating network and seeds channels based on the combination of all seed and network files provided")

        ch_seeds = Channel
            .fromPath(params.seeds.split(',').flatten(), checkIfExists: true)
            .combine(ch_network.map{meta, network -> meta.network_id})
            .map{seeds, network_id ->
                [ [ id: seeds.baseName + "." + network_id, seeds_id: seeds.baseName, network_id: network_id ] , seeds ]
            }

        // Add sp files, if provided (currently does not check if the number of the shortest paths matches the number of the networks and does not work with missing values)
        if(shortest_paths_param_set){
            ch_network = ch_network.merge(
                Channel
                .fromPath(params.shortest_paths.split(',').flatten())
            )
        } else{
            ch_network = ch_network.map{meta, network -> [meta, network, file("${projectDir}/assets/NO_FILE", checkIfExists: true)]}
        }

        // Add permuted network folders, if provided (currently does not check if the number of the shortest paths matches the number of the networks and does not work with missing values)
        if(permuted_networks_param_set){
            ch_network = ch_network.merge(
                Channel
                .fromPath(params.permuted_networks.split(',').flatten())
            )
        } else{
            ch_network = ch_network.map{meta, network, sp -> [meta, network, sp, []]}
        }

    } else {
        error("You need to specify either a sample sheet (--input) or the seeds (--seeds) and network (--network) files")
    }

    // check if IDs are unique
    ch_network.map{ meta, network, sp, permuted_networks -> [meta.id] }
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

    // separate network channel into network, shoretes_paths, and permuted_networks
    ch_shortest_paths = ch_network.map{meta, network, sp, permuted_networks ->
        [meta, sp.size() > 0 ? sp : file("${projectDir}/assets/NO_FILE", checkIfExists: true)]
    }

    ch_permuted_networks = ch_network.map{meta, network, sp, permuted_networks ->
        [meta, permuted_networks.size() > 0 ? file(permuted_networks+"/*.gt") : []]
    }

    ch_network = ch_network.map{meta, network, sp, permuted_networks -> [meta, network]}

    emit:
    versions    = ch_versions
    seeds       = ch_seeds                      // channel: [ val(meta[id,seeds_id,network_id]), path(seeds) ]
    network     = ch_network                    // channel: [ val(meta[id,network_id]), path(network) ]
    shortest_paths = ch_shortest_paths          // channel: [ val(meta[id,network_id]), path(shortest_paths) ]
    permuted_networks = ch_permuted_networks    // channel: [ val(meta[id,network_id]), [path(permuted_network)] ]
}

/*
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    SUBWORKFLOW FOR PIPELINE COMPLETION
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

def module_empty = [:]
def visualization_skipped = [:]
def drugstone_skipped = [:]
workflow PIPELINE_COMPLETION {

    take:
    email           //  string: email address
    email_on_fail   //  string: email address sent on pipeline failure
    plaintext_email // boolean: Send plain-text email instead of HTML
    outdir          //    path: Path to output directory where results will be published
    monochrome_logs // boolean: Disable ANSI colour codes in log output
    hook_url        //  string: hook URL for notifications
    multiqc_report  //  string: Path to MultiQC report
    module_empty_status          //  map: Empty/not empty status per module
    visualization_skipped_status //  map: Skipped/not skipped status per module
    drugstone_skipped_status     //  map: Skipped/not skipped status per module


    main:
    summary_params = paramsSummaryMap(workflow, parameters_schema: "nextflow_schema.json")
    def multiqc_reports = multiqc_report.toList()

    module_empty_status
        .map{
            id, status -> module_empty[id] = status
        }

    visualization_skipped_status
        .map{
            id, status -> visualization_skipped[id] = status
        }

    drugstone_skipped_status
        .map{
            id, status -> drugstone_skipped[id] = status
        }


    //
    // Completion email and summary
    //
    workflow.onComplete {
        if (email || email_on_fail) {
            completionEmail(
                summary_params,
                email,
                email_on_fail,
                plaintext_email,
                outdir,
                monochrome_logs,
                multiqc_reports.getVal(),
            )
        }

        logWarnings(monochrome_logs=monochrome_logs, module_empty=module_empty, visualization_skipped=visualization_skipped, drugstone_skipped=drugstone_skipped)
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
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    FUNCTIONS
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
*/

prepared_networks_url = "https://zenodo.org/records/15049754/files/"
network_map = [
    string_min900: "string.human_links_v12_0_min900",
    string_min700: "string.human_links_v12_0_min700",
    string_physical_min900: "string.human_physical_links_v12_0_min900",
    string_physical_min700: "string.human_physical_links_v12_0_min700",
    biogrid: "biogrid.4_4_242_homo_sapiens",
    hippie_high_confidence: "hippie.v2_3_high_confidence",
    hippie_medium_confidence:"hippie.v2_3_medium_confidence",
    iid: "iid.human",
    nedrex: "nedrex.reviewed_proteins_exp",
    nedrex_high_confidence: "nedrex.reviewed_proteins_exp_high_confidence",
]
id_space_map = [
    entrez: "Entrez",
    ensembl: "Ensembl",
    symbol: "Symbol",
    uniprot: "UniProtKB-AC",
]

//
// Check if the network is a prepared network or a file
//
def mapPreparedNetwork(network, id_space) {
    if (network_map.containsKey(network)) {
        return file("${prepared_networks_url}${network_map[network]}.${id_space_map[id_space]}.gt", checkIfExists: true)
    } else {
        return file(network, checkIfExists: true)
    }
}

//
// Read a tsv file and return its content as a list of groovy maps
//
def List<Map<String, String>> readTsvAsListOfMaps(file) {
    def lines = file.readLines()

    if (lines.size() < 2) {
        throw new IllegalArgumentException("TSV must have at least one header and one data row")
    }

    def headers = lines[0].split("\t")
    def result = []

    lines.tail().each { line ->
        def values = line.split("\t")
        if (values.size() != headers.size()) {
            throw new IllegalArgumentException("Mismatch between header and data line: $line")
        }

        def rowMap = [:]
        headers.eachWithIndex { header, idx ->
            rowMap[header] = values[idx]
        }
        result << rowMap
    }

    return result
}

//
// Create MultiQC tsv custom content from a list of values
//
def multiqcTsvFromList(tsv_data, header) {
    def tsv_string = ""
    if (tsv_data.size() > 0) {
        tsv_string += "${header.join('\t')}\n"
        tsv_string += tsv_data.join('\n')
    }
    return tsv_string
}


//
// Validate channels from input samplesheet
//
def validateInputSamplesheet(input) {
    def (metas, fastqs) = input[1..2]

    // Check that multiple runs of the same sample are of the same datatype i.e. single-end / paired-end
    def endedness_ok = metas.collect{ meta -> meta.single_end }.unique().size == 1
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
    // Convert  to a named map so can be used as with familiar NXF ${workflow} variable syntax in the MultiQC YML file
    def meta = [:]
    meta.workflow = workflow.toMap()
    meta["manifest_map"] = workflow.manifest.toMap()

    // Pipeline DOI
    if (meta.manifest_map.doi) {
        // Using a loop to handle multiple DOIs
        // Removing `https://doi.org/` to handle pipelines using DOIs vs DOI resolvers
        // Removing ` ` since the manifest.doi is a string and not a proper list
        def temp_doi_ref = ""
        def manifest_doi = meta.manifest_map.doi.tokenize(",")
        manifest_doi.each { doi_ref ->
            temp_doi_ref += "(doi: <a href=\'https://doi.org/${doi_ref.replace("https://doi.org/", "").replace(" ", "")}\'>${doi_ref.replace("https://doi.org/", "").replace(" ", "")}</a>), "
        }
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

def logWarnings(monochrome_logs=true, module_empty=[:], visualization_skipped=[:], drugstone_skipped=[:]) {
    def colors = logColours(monochrome_logs)

    def module_empty_count = module_empty.count  { key, value -> value == true }
    def visualization_skipped_count = visualization_skipped.count  { key, value -> value == true }
    def drugstone_skipped_count = drugstone_skipped.count  { key, value -> value == true }

    if (workflow.success) {
        if (module_empty_count > 0) {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} Please check MultiQC report: ${module_empty_count}/${module_empty.size()} modules were empty and not used in subsequent processes.${colors.reset}-"
        }
        if (visualization_skipped_count > 0) {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} Please check MultiQC report: ${visualization_skipped_count}/${visualization_skipped.size()} modules were too large for visualization (> ${params.visualization_max_nodes} nodes). You can adjust the threshold with the '--visualization_max_nodes' parameter.${colors.reset}-"
        }
        if (drugstone_skipped_count > 0) {
            log.info "-${colors.purple}[$workflow.manifest.name]${colors.red} Please check MultiQC report: ${drugstone_skipped_count}/${drugstone_skipped.size()} modules were too large for Drugst.One (> ${params.drugstone_max_nodes} nodes). You can adjust the threshold with the '--drugstone_max_nodes' parameter.${colors.reset}-"
        }
    }
}

