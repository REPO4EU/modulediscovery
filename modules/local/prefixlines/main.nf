
process PREFIXLINES {
    tag "$meta.id"
    label 'process_single'

    input:
    tuple val(meta), path(file)
    val prefix


    output:
    tuple val(meta), path("${file.baseName}.prefixed.${file.extension}")

    when:
    task.ext.when == null || task.ext.when


    script:
    """
    sed -e 's/^/${prefix}/' $file > ${file.baseName}.prefixed.${file.extension}
    """
}
