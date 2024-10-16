#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir = 'atavide.out'

params.dbUrl = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
params.temp = "./tmp"

process downloadDatabases {
    publishDir "${params.outdir}/databases", mode: 'copy', mkdirs: true

    input:
    val dbUrl

    output:
    path 'GCF_000001405.40_GRCh38.p14_genomic.fna.gz'

    script:
    """
    wget -O GCF_000001405.40_GRCh38.p14_genomic.fna.gz ${dbUrl}
    """
}

process UniRef50DB{
    publishDir "${params.outdir}/databases", mode: 'copy', mkdirs: true

    output:
    path 'UniRef50'

    conda 'bioconda::mmseqs2=14.7e284'

    script:
    """
    mmseqs databases --threads ${task.cpus} UniRef50 UniRef50 ${params.temp}
    """
}

workflow install_workflow {
    downloadDatabases(params.dbUrl)
    UniRef50DB()
}

// Main workflow entry point
workflow {
    install_workflow()
}