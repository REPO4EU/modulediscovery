
include { BIOPAX_PARSER         } from '../../../modules/local/biopax/parser/main'
include { BIOPAX_VALIDATOR      } from '../../../modules/local/biopax/validator/main'

workflow GT_BIOPAX {

    take:                                   // Workflow inputs
    ch_modules
    idspace

    main:

    ch_versions = Channel.empty()                                           // For collecting tool versions

    BIOPAX_PARSER(ch_modules, idspace)
    ch_versions = ch_versions.mix(BIOPAX_PARSER.out.versions)

    BIOPAX_VALIDATOR(BIOPAX_PARSER.out.biopax.collect())
    ch_versions = ch_versions.mix(BIOPAX_VALIDATOR.out.versions)


    emit:
    module   = BIOPAX_VALIDATOR.out.validation       // channel: [ module ]              emit the module discovered by DIAMOnD
    versions = ch_versions              // channel: [ versions.yml ]        emit collected versions
}
