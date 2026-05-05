# kentsislab/proteomegenerator3: Changelog

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/)
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [1.3.0] - 2026-05-05

### Added

- LRAA (Long Read Assembly and Annotation) as an alternative long-read assembler via `--long_read_assembler lraa`
- StringTie long-read mode as an alternative long-read assembler via `--long_read_assembler stringtie`
- CLI parameter `--long_read_assembler` to select the long-read transcript assembler (options: `bambu`, `lraa`, `stringtie`)
- CLI parameter `--skip_lraa_discovery` to skip LRAA assembly and use pre-computed GTFs from the samplesheet (proceeds directly to merge/reannotate/quantify)
- Strand-specific ORF prediction with Transdecoder as the default
- FAI index support for reference genome

### Changed

- CI profile ordering (`docker,test` → ensures test resource limits take precedence on GitHub runners)
- Updated LRAA to version 0.16.1

### Fixed

- Resource limits not being applied on GitHub runners due to profile ordering
- Output naming for non-NDR assemblers
- Glob patterns on LRAA modules
- Duplicate transcript structures and pruning of transcripts exceeding chromosome boundaries
- Novel transcript annotation in reannotate GTF
- Don't re-run quantification on reannotated LRAA GTF
- Don't run Bambu in preprocessing if not selected as the assembler

## [1.2.2] - 2026-04-16

### Added

- CLI parameter `--min_orf_len` to modulate minimum ORF length for Transdecoder (default: 100)

### Fixed

- Fixed SEQKIT_RMDUP publish renaming bug for multiple fasta outputs; uses `ext.prefix` to set output filename instead of `saveAs` rename
- Fixed publish results when only long-read data is supplied (no short-reads)
- Fixed merging long- and short-read data in single sample mode
- Fixed publishDir evaluation in SEMERGE at runtime
- Propagate NDR in file naming for correct output directory structure
- Give short/long-read merged proteomes subject ID prefix for proper sample identification

## [1.2.1] - 2026-03-04

### Fixed

- Fixed `meta.tool` not being set for bambu assemblies when `--short_reads` is disabled, which prevented TRANSDECODER2FASTA from publishing output correctly

## [1.2.0] - 2026-02-10

### Added

- Short-read RNA-seq support via StringTie integration
  - New `--short_reads` flag to enable short-read transcript assembly and quantification
  - StringTie-based assembly runs in parallel with Bambu for long-read data
  - Short-read and long-read ORF predictions are merged into unified proteome database
- Local seqkit_rmdup module with additional outputs for duplicate sequence tracking (`-D` and `-d` flags)
- New `--fusions` flag to control inclusion of fusion predictions in final proteome database
- Stub test implementations for all local modules (BAMBU_ASSEMBLY, BAMBU_FILTER, BAMBU_READCLASSES, SEMERGE, TRANSDECODER2FASTA, FUSIONFASTA, MERGEFUSIONS)
- Conditional fusion processing logic that handles samples without fusion data

### Changed

- **BREAKING**: Samplesheet format updated to support multiple samples per subject/patient
  - Added required `subject_id` column for patient/subject identifier
  - Renamed `sample` column to `sample_id`
  - `subject_id` is available in metadata as `meta.subject_id` for future subject-level processing
- FASTA_MERGE_ANNOTATE subworkflow now handles both bambu and stringtie outputs
  - Branching logic to separate outputs by tool type
  - CAT_CAT_SAMPLES process to merge bambu and stringtie fastas per sample
  - Updated publishDir configuration for proper routing of outputs to NDR directories
- FASTA_MERGE_ANNOTATE subworkflow now accepts `run_fusions` parameter to conditionally process fusion data
- Fusion processing is now opt-in via `--fusions` flag rather than automatic when fusion files are present
- Improved workflow logic to prevent empty channel errors when fusion data is not provided
- Moved seqkit/rmdup to local modules to support custom duplicate tracking outputs
- README updated with detailed `--fusions` flag usage examples and requirements

### Fixed

- Empty channel error in FASTA_MERGE_ANNOTATE when fusion files are not provided
- Pipeline now correctly handles mixed datasets where some samples have fusions and others don't
- Proper publishing of stringtie outputs to appropriate NDR directories (using params.NDR fallback)

### Removed

- Deprecated WRITEFASTA module (replaced by TRANSDECODER2FASTA)

## [1.1.1] - 2025-11-14

### Changed

- `ci.yml` updated with disk cleanup

### Fixed

- Fixed bug in transdecoder modules preventing correct processing of multiple samples simultaneously (moved transdecoder/longorf to local module and corrected input tuple structure in transdecoder/predict to properly pass folder parameter)

## [1.1.0] - 2025-11-12

### Added

- Support for ctat-lr-fusion input format (replacing JAFFAL for fusion analysis)
- New module: MERGEFUSIONS for processing and merging fusion contigs
- New module: FUSIONFASTA for extracting fusion sequences
- SwissProt concatenation functionality for multi-sample proteome databases
- SEQKIT modules (rmdup and stats) for sequence deduplication and statistics
- `--skip_multisample` flag to skip multi-sample transcript merging
- Support for starting from cached Bambu read class files via `rcFile` column in samplesheet
- BLAST integration into Transdecoder workflow for improved ORF prediction
- CLAUDE.md documentation for AI-assisted development guidance

### Changed

- Fusion workflow migrated from JAFFAL to ctat-lr-fusion format
- Enhanced transdecoder2fasta module with improved FASTA formatting
- Improved FASTA merge and annotation workflow
- Updated test configurations and snapshots

### Fixed

- Various fixes for predict_orfs subworkflow
- Improved samtools index handling
- Parameter and schema validation fixes
- PublishDir path corrections

## [1.0.0] - 2025-11-01

First stable release of kentsislab/proteomegenerator3.

### Changed

- Updated README.md with version 1.0.0 and biorxiv citation.

## [1.0.0dev] - 2025-07-23

Initial release of kentsislab/proteomegenerator3, created with the [nf-core](https://nf-co.re/) template.

### `Added`

- Transcript assembly & quant with bambu
- Flags to pre-filter reads before assembly on read length & mapq (useful for samples with inadequate QC)
- Flags to pre-filter reads on accessory chromosomes (can sometimes cause issues for Bambu)
- Flag to adjust NDR in bambu
- ORF prediction for transcripts & fusions using TransDecoder
- reformatting of fasta for use with MSFragger, DIA-NN, and Spectronaut
- nf-test and test datasets
- updated README.md
