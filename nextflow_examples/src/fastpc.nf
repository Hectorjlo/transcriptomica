// Define main parameters
params.fqs_dir = ""
params.qc_out = ""
params.fastp_dir = ""
params.cores = "2"

// Define the process to execute fastqc
process fastqc {
    // Output files of the process will be moved in "params.fastp_dir" 
    // due to file size
    publishDir params.qc_out, mode: 'move' 
    
    // Arguments expected by the process
    input:
    tuple val(SRR_name), path(fastq_path)
    val cores 
    val phase
    // Output expected by the process
    output:
    path "${SRR_name}_${phase}_fqc"

    // fastqc run
    script:
    """
    mkdir -p ${SRR_name}_${phase}_fqc
    /export/apps/bioconda/bin/fastqc ${fastq_path} -o ${SRR_name}_${phase}_fqc -t ${cores}
    """
}

process trim_by_fastp {
    // Output of the process will be moved due to file size
    publishDir params.fastp_dir, mode: 'move'

    // Arguments expected
    input:
    tuple val(SRR_name), path(fastq_path)
    val cores

    // Output expected
    output:
    tuple val(SRR_name), path("clean_*.fastq")

    // fastp run
    script:
    """
    /export/apps/bioconda/envs/fastp/bin/fastp -i ${fastq_path[0]} -I ${fastq_path[1]} \
          -o clean_${SRR_name}_1.fastq -O clean_${SRR_name}_2.fastq \
          --cut_front --cut_front_window_size 4 --cut_front_mean_quality 20 \
          --cut_tail --cut_tail_window_size 4 --cut_tail_mean_quality 20 \
          --detect_adapter_for_pe --trim_poly_g -l 50 -w ${cores}
    """
}

workflow raw_qc {
    take:
    SRR_pairs

    main:
    fastqc(SRR_pairs, params.cores, "raw")
}

workflow clean_qc {
    take:
    SRR_pairs

    main:
    fastqc(SRR_pairs, params.cores, "cleaned")
}

workflow {
    main:
    SRR_pairs = channel.fromFilePairs("${params.fqs_dir}/*_{1,2}.fastq")

    raw_qc(SRR_pairs)
    trimmed = trim_by_fastp(SRR_pairs, params.cores)

    clean_qc(trimmed)
}