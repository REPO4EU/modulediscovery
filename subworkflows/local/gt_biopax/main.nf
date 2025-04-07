
include { BIOPAX_PARSER         } from '../../../modules/local/biopax/parser/main'
include { BIOPAX_VALIDATOR      } from '../../../modules/local/biopax/validator/main'

workflow GT_BIOPAX {

    take:                                   // Workflow inputs
    ch_modules
    idspace
    validate_online

    main:

    ch_versions = Channel.empty()                                           // For collecting tool versions

    BIOPAX_PARSER(ch_modules, idspace, params.add_variants)                        // Parse the biopax files
    ch_versions = ch_versions.mix(BIOPAX_PARSER.out.versions)

    if (!params.add_variants){
        BIOPAX_VALIDATOR(BIOPAX_PARSER.out.biopax.collect(), validate_online) // Validate the biopax files
        ch_versions = ch_versions.mix(BIOPAX_VALIDATOR.out.versions)
        module = BIOPAX_VALIDATOR.out.validation
    } else {
        module = Channel.empty()
    }


    emit:
    module   = module
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
