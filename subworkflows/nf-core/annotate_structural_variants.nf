//
// A subworkflow to annotate structural variants.
//

include { SVDB_QUERY     } from '../../modules/nf-core/modules/svdb/query/main'
include { PICARD_SORTVCF } from '../../modules/nf-core/modules/picard/sortvcf/main'

workflow ANNOTATE_STRUCTURAL_VARIANTS {

    take:
        vcf         // channel: [ val(meta), path(vcf) ]
        sv_dbs      // file: dbs.csv
        fasta       // file: genome.fasta
        seq_dict    // file: genome.dict

    main:
        ch_versions = Channel.empty()

        Channel.fromPath(sv_dbs)
            .splitCsv ( header:true )
            .multiMap { row ->
                vcf_dbs:  row.filename
                in_frqs:  row.in_freq_info_key
                in_occs:  row.in_allele_count_info_key
                out_frqs: row.out_freq_info_key
                out_occs: row.out_allele_count_info_key
            }
            .set { ch_svdb_dbs }

        SVDB_QUERY(vcf,
            ch_svdb_dbs.in_occs.toList(),
            ch_svdb_dbs.in_frqs.toList(),
            ch_svdb_dbs.out_occs.toList(),
            ch_svdb_dbs.out_frqs.toList(),
            ch_svdb_dbs.vcf_dbs.toList()
            )

        PICARD_SORTVCF(SVDB_QUERY.out.vcf,
            fasta,
            seq_dict
        )

    emit:
        versions               = ch_versions.ifEmpty(null)      // channel: [ versions.yml ]
}