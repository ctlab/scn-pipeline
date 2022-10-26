# Changelog

## [1.0.4](https://github.com/ctlab/scn-pipeline/compare/v1.0.3...v1.0.4) (2022-10-26)


### Bug Fixes

* added priorities to run DAG depth-first ([605924e](https://github.com/ctlab/scn-pipeline/commit/605924e8a3df99a7a1e79db302e9429f6640d2c4))
* now ignoring bam index files, if they are present. fixing error where CB tag is not present in the first read of bam file ([cc886d3](https://github.com/ctlab/scn-pipeline/commit/cc886d3c18092733729d6d9a0f017bc660285e3d))
* now refering to factor variable explicitly, fixes [#3](https://github.com/ctlab/scn-pipeline/issues/3) ([29ee892](https://github.com/ctlab/scn-pipeline/commit/29ee892d698d544740b85b4ec805cc0f8c8e53f9))
* regular expression now works better with dots in the filename ([df86c70](https://github.com/ctlab/scn-pipeline/commit/df86c704b3b8d10c50513740289f24b712144b34))
* removed check for barcode length in bam files ([1972d1a](https://github.com/ctlab/scn-pipeline/commit/1972d1a117822b5a0501e7c222beacdcbce6a0c5))
* rename read_length to barcode_length, was causing inability to find reliably find index reads ([2258e89](https://github.com/ctlab/scn-pipeline/commit/2258e894747412fd6ac2d20ad60b6303fe1af713))
* set resource priority higher than data downloading, so it starts earlier ([53d99aa](https://github.com/ctlab/scn-pipeline/commit/53d99aa28c444060cc1bd1aac13ab459bc58604a))

## [1.0.3](https://github.com/ctlab/scn-pipeline/compare/v1.0.2...v1.0.3) (2022-10-04)

### Refactoring

* Now pipeline is compliant with `snakemake --lint`
* Many fixes to utilise relative paths ([edf2dec](https://github.com/ctlab/scn-pipeline/commit/edf2dec5e530a54dabdcd83ccaa1d780aed96461), [78c9b8b](https://github.com/ctlab/scn-pipeline/commit/78c9b8b8a1a36e5098eeddb4148f304901b43948), [7404967](https://github.com/ctlab/scn-pipeline/commit/740496761aca12019419b2d440276a837a97ddda), [0d25008](https://github.com/ctlab/scn-pipeline/commit/0d2500810baa458180e8e1d717360d9413912e54), and other commits)
 
### Bug Fixes
* Removed unused files - `bam_to_fastq` and `reports.smk` ([18b8a09](https://github.com/ctlab/scn-pipeline/commit/18b8a093bdea547389a654936e3b3c18451bbe98), [b16690a](https://github.com/ctlab/scn-pipeline/commit/b16690aa3e2b31679658ca780bebbde1e80033b7))
* Fixes to calculation of markers ([351b0ea](https://github.com/ctlab/scn-pipeline/commit/351b0ea53ec224ea809a73eeddd03245d3fc1deb))


## [1.0.2](https://github.com/ctlab/scn-pipeline/compare/v1.0.1...v1.0.2) (2022-09-20)


### Bug Fixes

* fixes to documentation to make it more transparent ([18c3e0e](https://github.com/ctlab/scn-pipeline/commit/18c3e0e0901bbdc43d196b151f8297db249b8fac))

## [1.0.1](https://github.com/ctlab/scn-pipeline/compare/v1.0.0...v1.0.1) (2022-09-02)


### Bug Fixes

* now ncbi_prefetch is not required for define_tech since, prefetch is now required for parallel-fastq-dump with small enough X, fixes [#2](https://github.com/ctlab/scn-pipeline/issues/2) ([fe24d74](https://github.com/ctlab/scn-pipeline/commit/fe24d74868e5521561992ebbb7c5f405a5a971c4))

## 1.0.0 (2022-09-01)

Releasing the scn-pipeline v.1.0.0

Currently supported technology - 10x and DropSeq


### Features

* added documentation. preparing for realease ([a30c03d](https://github.com/ctlab/scn-pipeline/commit/a30c03dc2d8f3c529899e75423c33af4f651f762))
* added release-please github action ([8081cda](https://github.com/ctlab/scn-pipeline/commit/8081cda222f16da2a55bf65f250ac320f80b1c04))
