#!/usr/bin/env nextflow

/*
 * eveHap - mtDNA Haplogroup Classification Pipeline
 *
 * Example Nextflow workflow for batch haplogroup classification.
 * Designed for HPC and cloud environments.
 *
 * Usage:
 *   nextflow run main.nf --input 'samples/*.bam' --tree rsrs --reference rcrs --outdir results/
 *   nextflow run main.nf --input 'samples/*.vcf.gz' --tree /path/to/tree.json --reference /path/to/ref.fasta
 */

nextflow.enable.dsl=2

// Parameters - tree and reference are REQUIRED
params.input = null                          // Input files (glob pattern)
params.outdir = 'results'                    // Output directory
params.tree = null                           // REQUIRED: 'rsrs', 'rcrs', or path to file
params.reference = null                      // REQUIRED: 'rcrs' or path to file
params.method = 'auto'                       // Classification method
params.damage_filter = false                 // Enable damage filtering
params.sample_id = null                      // Sample ID for multi-sample VCF

// Validate inputs
if (!params.input) {
    error """
    ════════════════════════════════════════════════════════════════════
    ERROR: --input is required

    Please provide input files:
      nextflow run main.nf --input 'samples/*.bam' --tree rsrs --reference rcrs
    ════════════════════════════════════════════════════════════════════
    """
}

if (!params.tree) {
    error """
    ════════════════════════════════════════════════════════════════════
    ERROR: --tree is required

    Options:
      --tree rsrs              Use bundled RSRS tree (~5400 haplogroups)
      --tree rcrs              Use bundled rCRS tree (~2400 haplogroups)
      --tree /path/to/tree.json   Use custom tree file

    To download the mitoLeaf tree (6400+ haplogroups):
      evehap download --outdir ./evehap_data/
      nextflow run main.nf --tree ./evehap_data/mitoleaf_tree.json ...
    ════════════════════════════════════════════════════════════════════
    """
}

if (!params.reference) {
    error """
    ════════════════════════════════════════════════════════════════════
    ERROR: --reference is required

    Options:
      --reference rcrs         Use bundled rCRS reference
      --reference /path/to/reference.fasta   Use custom reference

    To download references:
      evehap download --outdir ./evehap_data/
    ════════════════════════════════════════════════════════════════════
    """
}

// Log parameters
log.info """
╔═══════════════════════════════════════════════╗
║           eveHap Haplogroup Pipeline          ║
╚═══════════════════════════════════════════════╝

Input files  : ${params.input}
Output dir   : ${params.outdir}
Tree         : ${params.tree}
Reference    : ${params.reference}
Method       : ${params.method}
Damage filter: ${params.damage_filter}
"""

// Processes
process EVEHAP_CLASSIFY {
    tag "${sample_id}"
    publishDir "${params.outdir}", mode: 'copy'

    input:
    tuple val(sample_id), path(input_file)

    output:
    tuple val(sample_id), path("${sample_id}.haplogroup.tsv"), emit: results
    path "${sample_id}.haplogroup.json", emit: json

    script:
    def tree_opt = "--tree ${params.tree}"
    def ref_opt = "--reference ${params.reference}"
    def method_opt = "--method ${params.method}"
    def damage_opt = params.damage_filter ? "--damage-filter" : ""
    def sample_opt = params.sample_id ? "--sample-id ${params.sample_id}" : ""

    """
    evehap classify \\
        ${input_file} \\
        ${tree_opt} \\
        ${ref_opt} \\
        ${method_opt} \\
        ${damage_opt} \\
        ${sample_opt} \\
        --output-format tsv \\
        -o ${sample_id}.haplogroup.tsv \\
        -q

    evehap classify \\
        ${input_file} \\
        ${tree_opt} \\
        ${ref_opt} \\
        ${method_opt} \\
        ${damage_opt} \\
        ${sample_opt} \\
        --output-format json \\
        -o ${sample_id}.haplogroup.json \\
        -q
    """
}

process EVEHAP_DAMAGE {
    tag "${sample_id}"
    publishDir "${params.outdir}/damage", mode: 'copy'

    input:
    tuple val(sample_id), path(bam_file)

    output:
    tuple val(sample_id), path("${sample_id}.damage.txt"), emit: damage

    when:
    bam_file.extension == 'bam' || bam_file.extension == 'cram'

    script:
    """
    evehap damage ${bam_file} > ${sample_id}.damage.txt
    """
}

process MERGE_RESULTS {
    publishDir "${params.outdir}", mode: 'copy'

    input:
    path tsv_files

    output:
    path "all_haplogroups.tsv"

    script:
    """
    # Merge all TSV files, keeping header from first file only
    head -1 \$(ls *.haplogroup.tsv | head -1) > all_haplogroups.tsv
    for f in *.haplogroup.tsv; do
        tail -n +2 "\$f" >> all_haplogroups.tsv
    done
    """
}

// Workflow
workflow {
    // Create input channel
    Channel
        .fromPath(params.input)
        .map { file ->
            def sample_id = file.baseName.replaceAll(/\.(chrM|chrMT|MT)$/, '')
            tuple(sample_id, file)
        }
        .set { input_ch }

    // Run classification
    EVEHAP_CLASSIFY(input_ch)

    // Run damage analysis for BAM files
    input_ch
        .filter { sample_id, file ->
            file.extension == 'bam' || file.extension == 'cram'
        }
        .set { bam_ch }
    EVEHAP_DAMAGE(bam_ch)

    // Merge results
    EVEHAP_CLASSIFY.out.results
        .map { sample_id, tsv -> tsv }
        .collect()
        .set { all_tsvs }

    MERGE_RESULTS(all_tsvs)
}

// Completion handler
workflow.onComplete {
    log.info """
    ════════════════════════════════════════════════
    Pipeline completed!
    Results: ${params.outdir}/all_haplogroups.tsv
    ════════════════════════════════════════════════
    """.stripIndent()
}
