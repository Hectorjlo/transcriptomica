params.message = ""
params.outdir = ""

process underscore {
    input:
    val message

    output:
    path "exit.txt"

    script:
    """
    echo "$message" | tr " " "_" > exit.txt
    """
}

process addon {
    publishDir params.outdir, mode: 'copy'

    input:
    path tr_in

    output:
    path "movemo.txt"
    
    script:
    """
    cat "$tr_in" > movemo.txt
    """
}
workflow {
    raiz = channel.of(params.message)
    brazo = underscore(raiz)
    addon(brazo)
}