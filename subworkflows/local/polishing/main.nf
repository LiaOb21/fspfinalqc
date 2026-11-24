include { MINIBWA_INDEX } from '../../../modules/nf-core/minibwa/index/main'
include { MINIBWA_MAP   } from '../../../modules/nf-core/minibwa/map/main'
include {SAMTOOLS_VIEW  } from '../../../modules/nf-core/samtools/view/main'
include {SAMTOOLS_SORT  } from '../../../modules/nf-core/samtools/sort/main'
include {SAMTOOLS_FASTQ } from '../../../modules/nf-core/samtools/fastq/main'


workflow POLISHING {

    take:
    ch_samplesheet // channel samplesheet: [[id:sample, genus, family, order, class, phylum, single_end:false], [sample_R1.fastp.fastq.gz, sample_R2.fastp.fastq.gz]]
    ch_samplesheet_fasta // channel fasta: [[id:sample, genus, family, order, class, phylum, single_end:false], sample_corrected.fasta]
    ch_samplesheet_reads_fasta // channel reads_fasta: [[id:sample, genus, family, order, class, phylum, single_end:false], [sample_R1.fastp.fastq.gz, sample_R2.fastp.fastq.gz], sample_corrected.fasta]

    main:
    // minibwa index indexes the fasta file and has as input tuple val(meta), path(fasta). We can directly use ch_samplesheet_fasta
    MINIBWA_INDEX ( ch_samplesheet_fasta )
    // minibwa map maps the reads to the reference and has as input:
    // tuple val(meta), path(reads)
    // tuple val(meta2), path(index)
    // tuple val(meta3), path(fasta)

    ch_map_input = ch_samplesheet_reads_fasta
        .join( MINIBWA_INDEX.out.index, failOnMismatch: true )   // [ meta, reads, fasta, index ]
        .multiMap { meta, reads, fasta, index ->
            reads: [ meta, reads ]
            index: [ meta, index ]
            fasta: [ meta, fasta ]
        }

    MINIBWA_MAP ( ch_map_input.reads, ch_map_input.index, ch_map_input.fasta )

    SAMTOOLS_SORT ( ch_bam )

    SAMTOOLS_INDEX ( SAMTOOLS_SORT.out.bam )

    emit:
    // TODO nf-core: edit emitted channels
    bam      = SAMTOOLS_SORT.out.bam           // channel: [ val(meta), [ bam ] ]
    bai      = SAMTOOLS_INDEX.out.bai          // channel: [ val(meta), [ bai ] ]
    csi      = SAMTOOLS_INDEX.out.csi          // channel: [ val(meta), [ csi ] ]
}
