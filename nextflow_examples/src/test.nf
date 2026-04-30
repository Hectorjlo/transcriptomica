params.message = "og message"

process send_msg {
    input:
    val message

    script:
    """
    echo "nextflow dice: $message"
    """
}

workflow {
    main:
    datos = channel.of(params.message)
    send_msg(datos)
}