#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input = './testData/paired/'
params.outdir = 'atavide.out'

//databases 
params.host = './testData/host_reads.fna.gz'
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
        hostFilter(ch_trimmed, ch_genome)
    }
}
