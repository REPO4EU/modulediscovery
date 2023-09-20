
process PREFIXLINES {
    tag "$file"
    label 'process_single'

    input:
    path file
    val prefix


    output:
    path "${file.baseName}.prefixed.${file.extension}"

    when:
    task.ext.when == null || task.ext.when


    script:
    """
    sed -e 's/^/${prefix}/' $file > ${file.baseName}.prefixed.${file.extension}
    """
}
