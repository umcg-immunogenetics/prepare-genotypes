nextflow.enable.dsl=2

Channel
    .of(params.chromosomes)
    .set { chr_ch }

workflow {

    /*
     * Step 1: finalreport → lgen/map/fam
     */
    lgen_ch = CREATE_LGEN(
        file(params.input.final_report),
        file(params.input.sample_info)
    )

    /*
     * Step 2: lgen → plink bed
     */
    plink_ch = LGEN_TO_PLINK(lgen_ch)

    /*
     * Step 3: update build
     */
    updated_ch = UPDATE_BUILD(plink_ch)

    /*
     * Step 4: plink freq
     */
    freq_ch = PLINK_FREQ(updated_ch)

    /*
     * Step 5: HRC formatting
     */
    hrc_ch = MAKE_HRC_VCF(freq_ch)

    /*
     * Step 6: checkVCF (parallel per chr)
     */
    CHECK_VCF(chr_ch, hrc_ch)

    /*
     * Step 7: bgzip (parallel per chr)
     */
    BGZIP_VCF(chr_ch, hrc_ch)
}

process CREATE_LGEN {

    publishDir "${params.output_dir}", mode: 'copy'

    input:
    path final_report
    path sample_info

    output:
    tuple val(params.dataset),
          path("${params.dataset}.lgen"),
          path("${params.dataset}.map"),
          path("${params.dataset}.fam")

    script:
    """
    module load RPlus/4.2.1-foss-2022a-v22.12.1

    Rscript ${projectDir}/scripts/create_lgen.R \
        $final_report \
        $sample_info \
        ${params.dataset}.lgen \
        ${params.dataset}.map \
        ${params.dataset}.fam
    """
}

process LGEN_TO_PLINK {

    memory '16 GB'
    publishDir "${params.output_dir}", mode: 'copy'

    input:
    tuple val(prefix),
          path(lgen),
          path(map),
          path(fam)

    output:
    tuple val(prefix),
          path("${prefix}.bed"),
          path("${prefix}.bim"),
          path("${prefix}.fam")

    script:
    """
    module load PLINK/1.9-beta6-20190617

    plink --lfile ${prefix}
          --make-bed \
          --out ${prefix}
          --memory 16000
    """
}

process UPDATE_BUILD {

    process UPDATE_BUILD {

    publishDir "${params.output_dir}/GSA_updated", mode: 'copy'

    input:
    tuple val(prefix),
          path(bed),
          path(bim),
          path(fam)

    output:
    tuple val(prefix),
          path("${prefix}_GSA_updated.bed"),
          path("${prefix}_GSA_updated.bim"),
          path("${prefix}_GSA_updated.fam")

    script:
    """
    module load PLINK/1.9-beta6-20190617

    ${projectDir}/scripts/update_build.sh \
        ${prefix} \
        /groups/umcg-immunogenetics/tmp02/users/NathanRibeiro/tools/GSA_strand_v3/GSAMD-24v3-0-EA_20034606_A1-b37.strand \
        ${prefix}_GSA_updated
    """
}
}

process PLINK_FREQ {

    publishDir "${params.output_dir}/GSA_updated", mode: 'copy'

    input:
    path updated_prefix

    output:
    path "${params.dataset}_GSA_updated.frq*"

    script:
    """
    module load PLINK/1.9-beta6-20190617

    plink --freq \
          --bfile ${updated_prefix.baseName} \
          --out ${params.dataset}_GSA_updated
    """
}

process MAKE_HRC_VCF {

    publishDir "${params.output_dir}/HRC_formatted", mode: 'copy', pattern: "*-updated-chr*.vcf"

    input:
    path freq_prefix

    output:
    path "${params.dataset}_GSA_updated-updated-chr*.vcf", emit: vcf

    script:
    """
    module load foss/2022a
    module load PerlPlus/5.34.1-GCCcore-11.3.0-v22.11.1
    module load PLINK/1.9-beta6-20190617

    perl ${projectDir}/scripts/HRC-1000G-check-bim-NoReadKey.pl \
        -b ${freq_prefix.parent}/${params.dataset}_GSA_updated.bim \
        -f ${freq_prefix.parent}/${params.dataset}_GSA_updated.frq \
        -r ${params.HRC_reference} \
        -h

    sh Run-plink.sh
    """
}

process CHECK_VCF {

    publishDir "${params.output_dir}/checkVCF", mode: 'copy'

    input:
    val chr
    path hrc_vcf

    output:
    path "check-chr${chr}*"

    script:
    """
    export LD_LIBRARY_PATH=${params.openssl_path}:\$LD_LIBRARY_PATH

    ${params.python2_path}/python2.7 scripts/checkVCF.py \
        -r ${params.checkvcf_ref} \
        -o check-chr${chr} \
        ${hrc_vcf}
    """
}

process BGZIP_VCF {

    publishDir "${params.output_dir}/VCF_compressed", mode: 'copy'

    input:
    val chr
    path hrc_vcf

    output:
    path "${params.dataset}_final-chr${chr}.vcf.gz"

    script:
    """
    module load BCFtools/1.22-GCCcore-13.3.0

    bcftools view \
        ${hrc_vcf} \
        -Oz \
        -o ${params.dataset}_final-chr${chr}.vcf.gz
    """
}