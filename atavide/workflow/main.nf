#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input = './testData/paired/'
params.outdir = 'atavide.out'

//databases 
params.db = './atavide.out/databases/'
params.host = '/atavide.out/databases/GCF_000001405.40_GRCh38.p14_genomic.fna.gz'
params.uniref = '/home/edwa0468/Databases/mmseqs/UniRef50/'
params.temp = '/scratch/user/nala0006/tmp'

// Check if the database directory exists
if (!file(params.db).isDirectory() || file(params.db).list().length == 0) {
    include 'install.nf'
    install_workflow()
}

// Define the QC process
process FASTP {
    // define output directory
    publishDir "${params.outdir}/QC", pattern: "*_trimmed_*.fq.gz", mode: 'copy', mkdirs: true
    
    // input files
    input:
    tuple val(id), path(r1), path(r2)

    output:
    tuple val(id), path("${id}_trimmed_R1.fq.gz"), path("${id}_trimmed_R2.fq.gz") // Using tuple to pack R1 and R2 outputs

    conda 'bioconda::fastp=0.23.4'

    script:
    """
    # run pre-mapping QC
    fastp -i ${r1} -I ${r2} \
        -o "${id}_trimmed_R1.fq.gz" \
        -O "${id}_trimmed_R2.fq.gz"
    """
}

process hostFilter {
    publishDir "${params.outdir}/QC", pattern: "*host_removed*.fq.gz", mode: 'copy', mkdirs: true

    // input files
    input:
    tuple val(id), path(trimmed_r1), path(trimmed_r2) // Tuple to unpack the id, trimmed R1 and R2
    path genome // the genome input

    output:
    path "${id}_host_removed_R1.fq.gz"
    path "${id}_host_removed_R2.fq.gz"
    path "${id}_host_removed_singletons.fq.gz"

    conda 'bioconda::minimap2=2.28 bioconda::samtools=1.21'
 
    script:
    """
    minimap2 -t 32 --split-prefix=${params.temp} -a -xsr ${genome} ${trimmed_r1} ${trimmed_r2} | samtools view -bh - > "${params.temp}/output.bam"
    samtools index "${params.temp}/output.bam"
    samtools fastq -F 3584 -f 77 "${params.temp}/output.bam"  | gzip -c > "${id}_host_removed_R1.fq.gz"
    samtools fastq -F 3584 -f 141 "${params.temp}/output.bam"  | gzip -c > "${id}_host_removed_R2.fq.gz"
    samtools fastq -f 4 -F 1 "${params.temp}/output.bam"  | gzip -c > "${id}_host_removed_singletons.fq.gz"
    """
}

process readsTaxa {
    publishDir "${params.outdir}/mmseqs_taxonomy", mode: 'copy', mkdirs: true

    // input files
    input:
    tuple val(id), path(r1), path(r2) // Tuple to unpack the id, trimmed R1 and R2
    path uniref50db // mmseqs database

    output:
    path "${id}_lca.tsv"
    path "${id}_report"
    path "${id}_tophit_aln"
    path "${id}_tophit_report"

    conda 'bioconda::seqkit=2.8.2 bioconda::mmseqs2=14.7e284'
 
    script:
    """
    seqkit fq2fa ${r1} -o ${id}_R1.fasta
    seqkit fq2fa ${r2} -o ${id}_R2.fasta 
    mmseqs easy-taxonomy ${id}_R1.fasta ${id}_R2.fasta ${uniref50db}/UniRef50 ${id} ${params.temp} --threads 32
    """
}

process taxonkitLineage {
    publishDir "${params.outdir}/mmseqs_taxonomy", pattern: "*_taxonomy.tsv.gz", mode: 'copy', mkdirs: true

    input:
    tuple val(id), path(lca)

    output:
    path "${id}_taxonomy.tsv.gz"

    conda 'envs/taxonkit.yaml'

    script:
    """
    zcat ${lca} | awk '!s[\$2]++ {print \$2}' | taxonkit lineage | taxonkit reformat -P | awk -F"\t" -v OFS="\t" '{print \$1, \$3}' | gzip -c > ${id}_taxonomy.tsv.gz
    """
}

process add_taxonomy {
    publishDir "${params.outdir}/mmseqs_taxonomy", pattern: "*_lca_taxonomy.tsv.gz", mode: 'copy', mkdirs: true

    input:
    tuple val(id), path(tk), path(lca)

    output:
    path "${id}_lca_taxonomy.tsv.gz"

    conda 'envs/taxonkit.yaml'

    script:
    """
    python scripts/merge_taxonomy.py ${tk} ${lca} ${id}_lca_taxonomy.tsv.gz
    """
}

process summarise {
    publishDir "${params.outdir}/taxonomy", mode: 'copy', mkdirs: true

    input:
    tuple val(id), path(lcatax)

    output:
    path "${id}_summary.tsv.gz"

    script:
    """
    python scripts/summarise_taxonomy.py ${lcatax} ${id}_summary.tsv.gz
    """
}

process combine_all {
    publishDir "${params.outdir}/taxonomy_summary", mode: 'copy', mkdirs: true

    input:
    path summarised_files

    output:
    path "all_taxonomies.tsv"

    script:
    """
    cat ${summarised_files.join(' ')} > all_taxonomies.tsv
    """
}

workflow {
    // Channel for R1 reads
    ch_r1 = Channel.fromPath("${params.input}/*_R1*.fastq*")
        .map { file -> 
            def id = file.baseName.replace('_R1', '')
            return [id, file]
        }

    // Channel for R2 reads
    ch_r2 = Channel.fromPath("${params.input}/*_R2*.fastq*")
        .map { file -> 
            def id = file.baseName.replace('_R2', '')
            return [id, file]
        }

    // Join R1 and R2 channels by the sample id
    ch_paired = ch_r1.join(ch_r2, by: 0)

    // View for debugging (optional)
    ch_paired.view()

    // Run the FASTP process on paired files
    ch_trimmed = FASTP(ch_paired)

    // Channel for the UniRef50 database
    ch_uniref50 = Channel.fromPath(params.uniref)

    // Check if --hostFilter is passed
    if (params.hostFilter) {
        // Define the genome channel
        ch_genome = Channel.fromPath(params.host)

        // Run the hostFilter process with trimmed files and genome
        ch_filtered = hostFilter(ch_trimmed, ch_genome)

        // Set the input for readsTaxa to the output of hostFilter
        ch_taxa_output = readsTaxa(ch_filtered, ch_uniref50)
    } else {
        // If hostFilter is not run, use the outpqut from FASTP
        ch_taxa_output = readsTaxa(ch_trimmed, ch_uniref50)
    }

    
    ch_taxonomy = taxonkitLineage(ch_taxa_outputs.map ) // stuck here I need to take only one output from ch_taxa_output to this channel
    ch_addtaxa = add_taxonomy(ch_taxonomy)
    ch_summarise = summarise(ch_addtaxa)
    ch_summary = combine_all(ch_summarise)
}
