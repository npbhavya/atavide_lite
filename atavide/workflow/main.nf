#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input = './testData/paired/'
params.outdir = 'atavide.out'

//databases 
params.db = './atavide.out/databases'
params.host = './atavide.out/databases/GCF_000001405.40_GRCh38.p14_genomic.fna.gz'
params.temp = '/tmp/'

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
    minimap2 -t ${task.cpus} --split-prefix=${params.temp} -a -x sr ${genome} ${trimmed_r1} ${trimmed_r2} | samtools view -bh - > "${params.temp}/output.bam"
    samtools fastq -F 3584 -f 77 "${params.temp}/output.bam"  | gzip -c > "${id}_host_removed_R1.fq.gz"
    samtools fastq -F 3584 -f 141 "${params.temp}/output.bam"  | gzip -c > "${id}_host_removed_R2.fq.gz"
    samtools fastq -f 4 -F 1 "${params.temp}/output.bam"  | gzip -c > "${id}_host_removed_singletons.fq.gz"
    """
}

process readsTaxa {
    publishDir "${params.outdir}/mmseqs", pattern: "*_report.gz", mode: 'copy', mkdirs: true

    // input files
    input:
    tuple val(id), path(r1), path(r2) // Tuple to unpack the id, trimmed R1 and R2
    path UniRef50 // mmseqs database

    output:
    path "${id}_report.gz"
    path "${id}_lca.tsv.gz"
    path "${id}_tophit_report.gz"

    conda 'bioconda::mmseqs2=14.7e284'
 
    script:
    """
    atavide_lite/bin/fastq2fasta ${r1} ${id}_R1.fasta
    atavide_lite/bin/fastq2fasta ${r2} ${id}_R2.fasta
   
    # Check if FASTA files were created successfully
    if [ ! -f ${id}_R1.fasta ] || [ ! -f ${id}_R2.fasta ]; then
        echo "Error: FASTA files not created. Exiting."
        exit 1
    fi

    mmseqs easy-taxonomy ${id}_R1.fasta ${id}_R2.fasta ${UniRef50} ${id} ${params.temp} --start-sens 1 --sens-steps 3 -s 7 --threads ${task.cpus}
    """
}

workflow {
    // Channel for R1 reads
    ch_r1 = Channel.fromPath("${params.input}/*_R1*.fastq.gz")
        .map { file -> 
            def id = file.baseName.replace('_R1', '')
            return [id, file]
        }

    // Channel for R2 reads
    ch_r2 = Channel.fromPath("${params.input}/*_R2*.fastq.gz")
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

    // Check if --hostFilter is passed
    if (params.hostFilter) {
        // Define the genome channel
        ch_genome = Channel.fromPath(params.host)

        // Run the hostFilter process with trimmed files and genome
        ch_filtered = hostFilter(ch_trimmed, ch_genome)

        // Set the input for readsTaxa to the output of hostFilter
        ch_taxa_input = ch_filtered
    } else {
        // If hostFilter is not run, use the output from FASTP
        ch_taxa_input = ch_trimmed
    }
}
