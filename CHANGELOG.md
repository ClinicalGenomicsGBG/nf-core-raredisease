# nf-core/raredisease: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## v1.2.0 - Genie [xxxx-xx-xx]

### `Added`

- Use `nf-validation` plugin for parameter and samplesheet validation [#386](https://github.com/nf-core/raredisease/pull/386)

### `Changed`

- Tiddit updated to v3.6.1 [#385](https://github.com/nf-core/raredisease/pull/385)
- Installed the nf-core version of the sentieon/bwamemindex module [#397](https://github.com/nf-core/raredisease/pull/397)
- Installed the nf-core version of the sentieon/bwamem module [#398](https://github.com/nf-core/raredisease/pull/398)
- Installed the nf-core version of the sentieon/readwriter module [#399](https://github.com/nf-core/raredisease/pull/399)
- Installed the nf-core version of the sentieon/datametrics module [#400](https://github.com/nf-core/raredisease/pull/400)
- Installed the nf-core version of the sentieon/dedup module. The dedup module also contains a call to Sentieon's LocusCollector [#401](https://github.com/nf-core/raredisease/pull/401)
- Removing Sentieon-based BQSR. Recent Illumina sequencers tend to provide well-calibrated BQs, so BQSR may not provide much benefit [#402](https://github.com/nf-core/raredisease/pull/402)
- Installed the nf-core version of the sentieon/dnamodelapply module [#403](https://github.com/nf-core/raredisease/pull/403)
- Installed the nf-core version of the sentieon/wgsmetricsalgo module [#404](https://github.com/nf-core/raredisease/pull/404)
- Installed the nf-core version of the sentieon/dnascope module [#406](https://github.com/nf-core/raredisease/pull/406)

### `Fixed`

- Invalid GATK4 container which caused incorrect singularity downloads with nf-core download [nf-core/modules #3668](https://github.com/nf-core/modules/issues/3668)
- Make the default cram prefix same as markduplicates prefix [#392](https://github.com/nf-core/raredisease/pull/392)
- Sort ranked SV vcf before indexing with tabix [#393](https://github.com/nf-core/raredisease/pull/393)
- Make target bed file optional for WGS mode (Issue [#375](https://github.com/nf-core/raredisease/issues/375)) [#395](https://github.com/nf-core/raredisease/pull/395)
- Added constraints to block the pipeline from running CollectWgsMetrics on WES samples [#396](https://github.com/nf-core/raredisease/pull/396)

## v1.1.1 - Abu (Patch) [2023-07-26]

### `Fixed`

- Avoids errors thrown by bcftools concat due to sample names in input vcf files not being in same order [#388](https://github.com/nf-core/raredisease/pull/388)

## v1.1.0 - Abu [2023-07-21]

### `Added`

- Add GATK's cnv calling pipeline [#362](https://github.com/nf-core/raredisease/pull/362)
- GATK's ShiftFasta to generate all the files required for mitochondrial analysis [#354](https://github.com/nf-core/raredisease/pull/354)
- Feature to calculate CADD scores for indels [#325](https://github.com/nf-core/raredisease/pull/325)
- HmtNote to annotate mitochondria [#355](https://github.com/nf-core/raredisease/pull/355)
- MT del script to detect mitochondrial deletions [#349](https://github.com/nf-core/raredisease/pull/349)
- eKLIPse to identify large mitochondrial deletions [#365](https://github.com/nf-core/raredisease/pull/365)
- UPD+Chromograph to identify and visualize UPD sites and regions in the chromosomes [#364](https://github.com/nf-core/raredisease/pull/364) and [#366](https://github.com/nf-core/raredisease/pull/366)
- Added check for presence of case id for each sample in samplesheet [#357](https://github.com/nf-core/raredisease/pull/357)

### Fixed

- Avoiding publishing uncompressed VCF-file from `HMTNOTE_ANNOTATE`. (The corresponding compressed VCF-file still gets published.) [#368](https://github.com/nf-core/raredisease/pull/368)

## v1.0.0 - Aladdin [2023-06-01]

Initial release of nf-core/raredisease, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- FastQC read quality control
- Read mapping with BWAmem2/Sentieon
- Qualimap & Picard tools quality control metrics
- Call repeat expansions with ExpansionHunter and Stranger
- SNV calling with DeepVariant/Sentieon
- SV calling with Manta and TIDDIT
- SNV annotation with bcftools roh, vcfanno, and vep
- SV annotation with SVDB query and vep
- Separate workflow for analysing and annotating mitochondrial variants
- Call copy number variants in the SMN gene using SMNCopyNumberCaller
