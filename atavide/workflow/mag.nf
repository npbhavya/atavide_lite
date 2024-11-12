#!/usr/bin/env nextflow
nextflow.enable.dsl=2

params.input = 'atavide.out/QC/'
params.outdir = 'atavide.out'

params.gtdbtkDB= '/home/mall0133/miniconda3/envs/gtdbtk/share/gtdbtk-2.2.3/db'

// Define the QC process
process ASSEMBLY {
    // define output directory
    publishDir "${params.outdir}/assembly", mode: 'copy', mkdirs: true

    // input files
    input:
    tuple val(id), path(r1), path(r2)

    output:
    tuple val(id), path("${id}_spades/${id}_contigs.fasta")// Using tuple to pack R1 and R2 outputs

    conda 'bioconda::spades=4.0.0 bioconda::quast=5.2.0'

    script:
    """
    metaspades.py -1 ${r1} -2 ${r2} -o "${id}_spades" -t 32
    mv "${id}_spades"/contigs.fasta "${id}_spades"/"${id}"_contigs.fasta
    quast.py ${id}_spades/"${id}"_contigs.fasta -o ${id}_spades/quast
    """
}

process CONCATENATE {
    publishDir "${params.outdir}/vamb", mode: 'copy', mkdirs: true

    input:
    path (contigs)

    output:
    tuple path("allcontigs.mmi"), path("allcontigs.fa")

    conda 'bioconda::minimap2=2.28 bioconda::vamb=3.0.2 bioconda::samtools=1.21'

    script:
    def contigs_files = contigs.join(' ')
    """
    concatenate.py allcontigs.fa ${contigs_files}
    minimap2 -d allcontigs.mmi allcontigs.fa
    """
}

// Define the MINIMAP2 process
process MINIMAP2 {
    publishDir "${params.outdir}/vamb", mode: 'copy', mkdirs: true

    input:
    tuple val(id), path(r1), path(r2), path (index), path (fasta)

    output:
    path("${id}.bam")

    conda 'bioconda::minimap2=2.28 bioconda::vamb=3.0.2 bioconda::samtools=1.21'

    script:
    """
    minimap2 -t 32 -N 5 -ax sr ${index} --split-prefix mmsplit \\
        ${r1} ${r2} | samtools view -F 3584 -b --threads 32 > ${id}.bam
    """
}

process VAMB{
    publishDir "${params.outdir}/vamb", mode: 'copy', mkdirs: true

    input:
    tuple path (fasta), path (mmi), path (bams)

    output:
    path ("vamb-out/bins")

    conda 'bioconda::minimap2=2.28 bioconda::vamb=3.0.2 bioconda::samtools=1.21'

    script:
    // Collect all BAM files into a space-separated string
    """
    vamb --outdir vamb-out --fasta allcontigs.fa --bamfiles *.bam -o C --minfasta 200000 -t 32
    """
}

process GTDBTKANI {
    publishDir "${params.outdir}/vamb/gtdbtk", mode: 'copy', mkdirs:true

    input:
    path ("vamb-out/bins")

    output:
    path ("ani_rep")

    conda 'bioconda:: gtdbtk=2.4.0'

    script:
    """
    export GTDBTK_DATA_PATH=${params.gtdbtkDB}
    gtdbtk ani_rep --genome_dir vamb-out/bins -x fna --cpus 32
    """
}

process GTDBTK {
    publishDir "${params.outdir}/vamb/gtdbtk", mode: 'copy', mkdirs:true

    input:
    path ("vamb-out/bins")

    output:
    path ("classify")

    conda 'bioconda:: gtdbtk=2.4.0'

    script:
    """
    export GTDBTK_DATA_PATH=${params.gtdbtkDB}
    gtdbtk identify --genome_dir vamb-out/bins -x fna --cpus 32 --out_dir identify
    gtdbtk align --identify_dir identify --out_dir align --cpus 32
    gdbtk classify --genome_dir vamb-out/bins -x fna --cpus 32 --out_dir classify -f --align_dir align
    """
}

workflow {
    // Channel for R1 reads
    ch_r1 = Channel.fromPath("${params.input}/*_R1*")
        .map { file -> 
            def id = (file.baseName =~ /(.*)_R1/)[0][1]
            return [id, file]
        }

    // Channel for R2 reads
    ch_r2 = Channel.fromPath("${params.input}/*_R2*")
        .map { file -> 
            def id = (file.baseName =~ /(.*)_R2/)[0][1]
            return [id, file]
        }

    // Join R1 and R2 channels by the sample id
    ch_paired = ch_r1.join(ch_r2, by: 0)
    
    // Run assembly
    ch_contigs = ch_paired | ASSEMBLY | map { id, file -> file }
    
    // Collect contigs
    ch_combined = ch_contigs.collect()

    // Run concatenation
    ch_concat_result = CONCATENATE(ch_combined)
    
    // Run MINIMAP2 with paired reads and concatenated index
    ch_minimap_input = ch_paired.combine(ch_concat_result)

    //Run minimap2
    ch_minimap_output = MINIMAP2(ch_minimap_input)

    //Collect all the bamfiles
    ch_bam_list = ch_minimap_output.collect().map { bams -> 
        bams.unique() }

    // Prepare VAMB input - simplified version
    ch_vamb_input = ch_concat_result
        .combine(ch_bam_list.toList())
        .map {mmi, fasta, bams -> 
            tuple(fasta, mmi, bams)
        }

    // Run VAMB
    ch_vamb_output = VAMB(ch_vamb_input)

    //RUN GTDBTK
    ch_gtdb_ani = GTDBTKANI (ch_vamb_output)
    ch_gtdbtk = GTDBTK (ch_vamb_output)

}