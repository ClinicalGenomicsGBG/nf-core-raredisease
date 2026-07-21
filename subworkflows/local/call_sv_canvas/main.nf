include { BCFTOOLS_NORM           } from "../../../modules/nf-core/bcftools/norm/main"
include { CANVAS_GERMLINE         } from '../../../modules/nf-core/canvas/germline/main'
include { CANVAS_CREATESEG        } from '../../../modules/local/canvas_createseg/main'
include { GAWK as GAWK_CREATE_SEG } from '../../../modules/nf-core/gawk/main'
include { GAWK as GAWK_REFORMAT   } from '../../../modules/nf-core/gawk/main'
include { TABIX_BGZIPTABIX        } from "../../../modules/nf-core/tabix/bgziptabix/main"
include { TABIX_TABIX             } from "../../../modules/nf-core/tabix/tabix/main"

workflow CALL_SV_CANVAS {

    take:
    ch_bam_bai // channel: [meta, path(bam), path(bai) ]
    ch_snv     // channel: [meta, path(vcf)]
    ch_fasta   // channel: [meta, path(fasta)]
    ch_genomedir // channel: [meta, path(genomedir)]
    ch_m_ploidy_vcf // channel: [meta, path(male_ploidy_vcf)]
    ch_f_ploidy_vcf // channel: [meta, path(female_ploidy_vcf)]
    ch_filter_bed // channel: [meta, path(filter13)]
    reformat_vcf // boolean: [optional] [default: false] whether to reformat the canvas vcf output to match the expected nf-core/cnvkit vcf output

    main:

    ch_ploidy_vcf = ch_m_ploidy_vcf
      .combine(ch_f_ploidy_vcf)
      .map { meta1, male_vcf, _meta2, female_vcf -> tuple(meta1, male_vcf, female_vcf) }
      .first()

    CANVAS_GERMLINE(
        ch_bam_bai,
        ch_snv,
        ch_ploidy_vcf,
        ch_genomedir,
        ch_fasta,
        ch_filter_bed
    )

    GAWK_CREATE_SEG(
        CANVAS_GERMLINE.out.covandvarfreq,
        channel.value(file(moduleDir + '/bin/create_seg.awk'))
    )

    if (reformat_vcf) {
        BCFTOOLS_NORM(
            CANVAS_GERMLINE.out.vcf,
            ch_fasta
        )

        GAWK_REFORMAT (
            BCFTOOLS_NORM.out.vcf,
            channel.value(file(moduleDir + '/bin/reformat_canvas_vcf.awk'))
        )

        TABIX_BGZIPTABIX(
            GAWK_REFORMAT.out.output
        )
    } else {
        TABIX_TABIX(
            ch_vcf
        )

    }

    ch_vcf = reformat_vcf ? TABIX_BGZIPTABIX.out.gz_index.map 


    emit:
    vcf = ch_vcf
    index = TABIX_TABIX.out.index
    seg = GAWK_CREATE_SEG.out.output

}
