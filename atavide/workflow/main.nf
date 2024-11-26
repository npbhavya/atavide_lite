#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input = './testData/paired/'
params.outdir = 'atavide.out'

//databases 
params.db = './atavide.out/databases/'
params.host = './atavide.out/databases/GCF_000001405.40_GRCh38.p14_genomic.fna.gz'
params.uniref = '/home/edwa0468/Databases/mmseqs/UniRef50/'
params.temp = '/scratch/user/nala0006/tmp'
//params.taxonDB = "/home/edwa0468/Databases/NCBI/taxonomy/Oct_2024/"
params.taxonDB = "/home/nala0006/.taxonkit"

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
    fastp -i ${r1} -I ${r2} -o ${id}_trimmed_R1.fq.gz -O ${id}_trimmed_R2.fq.gz
    """
}

process HOSTFILTER {
    publishDir "${params.outdir}/QC", pattern: "*host_removed*.fq.gz", mode: 'copy', mkdirs: true

    // input files
    input:
    tuple val(id), path(trimmed_r1), path(trimmed_r2), path (genome) // the genome input

    output:
    tuple val(id), path ("${id}_host_removed_R1.fq.gz"), path ("${id}_host_removed_R2.fq.gz")

    conda 'bioconda::minimap2=2.28 bioconda::vamb=3.0.2 bioconda::samtools=1.21'
 
    script:
    """
    minimap2 -t 32 --split-prefix=${params.temp} -a -xsr ${genome} ${trimmed_r1} ${trimmed_r2} | samtools view -bh | samtools sort -o sorted_output.bam -
    #samtools sort output.bam -o sorted_output.bam
    samtools index sorted_output.bam
    samtools fastq -F 3584 -f 77 sorted_output.bam  | gzip -c > ${id}_host_removed_R1.fq.gz
    samtools fastq -F 3584 -f 141 sorted_output.bam  | gzip -c > ${id}_host_removed_R2.fq.gz
    samtools fastq -f 4 -F 1 sorted_output.bam  | gzip -c > ${id}_host_removed_singletons.fq.gz
    """
}

process MMSEQS {
    publishDir "${params.outdir}/mmseqs_taxonomy", mode: 'copy', mkdirs: true

    // input files
    input:
    tuple val(id), path(r1), path(r2), path (uniref50db) // mmseqs database

    output:
    tuple val(id), path ("${id}_lca.tsv"), path ("${id}_report"), path ("${id}_tophit_aln"), path ("${id}_tophit_report")

    conda 'bioconda::seqkit=2.8.2 bioconda::mmseqs2=14.7e284'
 
    script:
    """
    seqkit fq2fa ${r1} -o ${id}_R1.fasta
    seqkit fq2fa ${r2} -o ${id}_R2.fasta 
    mmseqs easy-taxonomy ${id}_R1.fasta ${id}_R2.fasta ${uniref50db}/UniRef50 ${id} ${params.temp} --threads 32
    """
}

process TAXONKIT {
    publishDir "${params.outdir}/mmseqs_taxonomy", pattern: "*_taxonomy.tsv.gz", mode: 'copy', mkdirs: true

    input:
    tuple val(id), path(lca), path (report), path (tophit_aln), path (tophit_report), path (taxonDB)

    output:
    tuple val(id), path(lca), path ("${id}_taxonomy.tsv.gz")

    conda 'bioconda::taxonkit=0.18.0'

    script:
    """
    export TAXONDB=${taxonDB}
    cat ${lca} | awk '!s[\$2]++ {print \$2}' | taxonkit lineage | taxonkit reformat -P | awk -F"\t" -v OFS="\t" '{print \$1, \$3}' | gzip -c > ${id}_taxonomy.tsv.gz
    """
}

process ADD_TAXONOMY {
    publishDir "${params.outdir}/mmseqs_taxonomy", pattern: "*_lca_taxonomy.tsv.gz", mode: 'copy', mkdirs: true

    input:
    tuple val(id), path(lca), path (tk), path (taxonDB)
    path script_path

    output:
    tuple val(id), path ("${id}_lca_taxonomy.tsv.gz")

    conda 'bioconda::taxonkit=0.18.0'

    script:
    """
    export tk=${tk}
    export lca=${lca} 
    export taxonDB=${taxonDB}
    export outputfile="${id}_lca_taxonomy.tsv.gz"
    python $script_path
    """
}

process SUMMARISE {
    publishDir "${params.outdir}/taxonomy", mode: 'copy', mkdirs: true

    input:
    tuple val(id), path(lcatax)
    path script_summary_path

    output:
    path ("${id}_summary")

    script:
    """
    export lcatax=${lcatax}
    export output="${id}_summary"
    python $script_summary_path
    """
}

// Define the COMBINE process
process COMBINE {
    // Specify where to save the output files
    publishDir "${params.outdir}/taxonomy", mode: 'copy', mkdirs: true

    input:
    tuple val (tax),  path(inputFiles)
    path script_join_path

    output:
    path "${tax}.raw.tsv.gz", emit: raw_output // Raw combined output file for the taxonomic rank
    path "${tax}.norm.tsv.gz", emit: norm_output // Normalized output file for the taxonomic rank

    script:
    """
    # Call a Python script to combine the input files into taxonomic-level summaries
    python ${script_join_path} ${inputFiles.collect { it + "/${tax}.tsv.gz" }.join(' ')} ${tax}.raw.tsv.gz ${tax}.norm.tsv.gz    
    """
}

process ANNOT {
    // Specify where to save the output files
    publishDir "${params.outdir}/annotation", mode: 'copy', mkdirs: true

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

    // Run the FASTP process on paired files
    ch_trimmed = FASTP(ch_paired)

    // Check if --hostFilter is passed
    if (params.hostFilter) {
        // Define the genome channel
        ch_genome = Channel.fromPath(params.host)
        ch_trimmed_input = ch_trimmed.combine(ch_genome)
        
        // Run the hostFilter process with trimmed files and genome
        ch_filtered = HOSTFILTER(ch_trimmed_input)

        //ch_filtered = HOSTFILTER(ch_trimmed, ch_genome)

        // Set the input for readsTaxa to the output of hostFilter
        ch_taxa_input = ch_filtered
    } else {
        // If hostFilter is not run, use the output from FASTP
        ch_taxa_input = ch_trimmed
    }
    
    // Channel for the UniRef50 database
    ch_uniref50 = Channel.fromPath(params.uniref)
    ch_taxa_inputWDB=ch_taxa_input.combine(ch_uniref50)
    ch_taxa_output = MMSEQS(ch_taxa_inputWDB)

    // TaxonDB set 
    ch_taxonDB = Channel.fromPath(params.taxonDB)
    ch_taxonkit_input=ch_taxa_output.combine(ch_taxonDB)
    ch_taxonomy = TAXONKIT(ch_taxonkit_input) 

    //set script to run 
    script_path = file("./atavide/scripts/merge_taxonomy.py")
    ch_addtaxa_input=ch_taxonomy.combine(ch_taxonDB)
    //ch_addtaxa_input.view()
    ch_addtaxa = ADD_TAXONOMY(ch_addtaxa_input, script_path)

    // set script to run 
    script_summary_path = file("./atavide/scripts/summarise_taxonomy.py")
    // Run the SUMMARISE process for each sample in ch_addtaxa
    ch_summarise = SUMMARISE(ch_addtaxa, script_summary_path).collect()

    script_join_path = file("./atavide/scripts/join.py")
    taxa = Channel.of ('phylum','kingdom','class','family','genus','species')
    ch_combine_input = taxa.combine(ch_summarise.toList())
    //ch_combine_input.view()
    ch_combine_output = COMBINE (ch_combine_input, script_join_path)
}
