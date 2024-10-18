#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.outdir = 'atavide.out'

params.dbUrl = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001/405/GCF_000001405.40_GRCh38.p14/GCF_000001405.40_GRCh38.p14_genomic.fna.gz"
// Define UniRef50 database URL for downloading
params.uniref50Url = 'https://ftp.uniprot.org/pub/databases/uniprot/uniref/uniref50/uniref50.fasta.gz'

params.temp = "/scratch/user/nala0006/tmp"

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

// Download UniRef50 database
process downloadUniRef50 {
    publishDir "${params.outdir}/databases", mode: 'copy', mkdirs: true

    output:
    path 'uniref50.fasta.gz'

    script:
    """
    # Download UniRef50 database
    wget ${params.uniref50Url} -O uniref50.fasta.gz
    """
}

process createUniRef50Db {
    publishDir "${params.outdir}/databases/", mode: 'copy', mkdirs: true

    input:
    path 'uniref50.fasta.gz'

    output:
    path 'UniRef50_db'

    conda 'bioconda::mmseqs2=14.7e284'

    script:
    """
    # Create database directory
    mkdir -p UniRef50_db
    
    # Create MMseqs2 database from UniRef50
    mmseqs createdb uniref50.fasta.gz UniRef50_db/UniRef50
    mmseqs createtaxdb UniRef50_db/UniRef50 UniRef50_db/taxonomy
    """
}

workflow install_workflow {
    downloadDatabases(params.dbUrl)
    downloadUniRef50()
    createUniRef50Db(downloadUniRef50.out)
}

// Main workflow entry point
workflow {
    install_workflow()
}