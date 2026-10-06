version = "2.1.1"
container_url = "xgrand/featurecounts:${version}"

params.fc_out = ""
params.fc_param = "2"      // strandness: 0 = unstranded, 1 = stranded, 2 = reversely stranded
params.fc_feature = "exon"
params.fc_attr = "gene_id"
params.fc_extra = ""       // args additionnels (e.g. "-p --countReadPairs")

process gff3_2_gtf {
    container = "dceoy/cufflinks"
    label "small_mem_mono_cpus"

    input:
    tuple val(genome_id), path(gff3_file)

    output:
    path "${genome_id}.gtf", emit: gtf

    script:
    """
    gffread ${gff3_file} -T -o ${genome_id}.gtf
    """
}

process featurecounts {
    container = "${container_url}"
    label "big_mem_multi_cpus"
    tag "$file_id"
    if (params.fc_out != "") {
        publishDir "results/${params.fc_out}", mode: 'copy'
    }

    input:
    tuple val(file_id), path(bam), path(bai)
    path(gtf)

    output:
    path "${file_id}.tsv", emit: counts

    script:
    def attr = params.fc_attr
    """
    featureCounts \\
        -T ${task.cpus} \\
        -a $gtf \\
        -F GTF \\
        -t ${params.fc_feature} \\
        -g ${attr} \\
        -s ${params.fc_param} \\
        ${params.fc_extra} \\
        -o ${file_id}.raw.tsv \\
        $bam

    # Comptages au format htseq-count : 2 colonnes, sans en-tête,
    # trié par gene_id, métriques QC (__*) en fin de fichier
    awk -F '\\t' 'NR > 2 { print \$1 "\\t" \$NF }' ${file_id}.raw.tsv \\
        | sort -k1,1 > ${file_id}.tsv

    # Métriques QC équivalentes aux pseudo-gènes htseq
    status_col=\$(awk -F '\\t' 'NR > 1 && \$1 == "Assigned" { print NF; exit }' ${file_id}.raw.tsv.summary)
    {
        printf "__no_feature\\t%s\\n" \$(awk -F '\\t' -v c=\$status_col '\$1 == "Unassigned_NoFeatures" { print \$c }' ${file_id}.raw.tsv.summary)
        printf "__ambiguous\\t%s\\n" \$(awk -F '\\t' -v c=\$status_col '\$1 == "Unassigned_Ambiguity" { print \$c }' ${file_id}.raw.tsv.summary)
        printf "__too_low_aQual\\t%s\\n" \$(awk -F '\\t' -v c=\$status_col '\$1 == "Unassigned_MappingQuality" { print \$c }' ${file_id}.raw.tsv.summary)
        printf "__not_aligned\\t%s\\n" \$(awk -F '\\t' -v c=\$status_col '\$1 == "Unassigned_Unmapped" { print \$c }' ${file_id}.raw.tsv.summary)
        printf "__alignment_not_unique\\t%s\\n" \$(awk -F '\\t' -v c=\$status_col '\$1 == "Unassigned_MultiMapping" { print \$c }' ${file_id}.raw.tsv.summary)
    } >> ${file_id}.tsv
    """
}

workflow featurecounts_with_gff {
    take:
    bam_tuple
    gff_file

    main:
    gff3_2_gtf(gff_file)
    featurecounts(bam_tuple, gff3_2_gtf.out.gtf)

    emit:
    counts = featurecounts.out.counts
}