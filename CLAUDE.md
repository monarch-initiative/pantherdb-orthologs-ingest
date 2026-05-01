# pantherdb-orthologs-ingest

This is a Koza ingest repository for transforming PantherDB ortholog data into Biolink model format.

## Project Structure

- `download.yaml` - Configuration for downloading PantherDB ortholog data
- `src/` - Transform code and configuration
  - `transform.py` / `transform.yaml` - Main transform for ortholog associations
  - `panther_orthologs_utils.py` - Utility functions for parsing gene info and mapping
  - `versions.py` - Per-ingest upstream version fetcher (consumed by `just metadata`)
- `scripts/` - Utility scripts (preprocessing, plus `write_metadata.py` which emits `output/release-metadata.yaml`)
- `tests/` - Unit tests for transforms
- `output/` - Generated nodes and edges (gitignored)
  - `release-metadata.yaml` - Per-build manifest of upstream sources, versions, artifacts (kozahub-metadata-schema)
- `data/` - Downloaded source data (gitignored)

## Key Commands

- `just run` - Full pipeline (download -> transform)
- `just download` - Download PantherDB ortholog data
- `just transform-all` - Run all transforms
- `just transform <name>` - Run specific transform
- `just metadata` - Emit `output/release-metadata.yaml`
- `just test` - Run tests

## Data Processing

This ingest processes the PantherDB AllOrthologs file, extracting orthology relationships between genes across multiple species. It uses an NCBI gene_info file to map species-specific identifiers to NCBIGene IDs when direct database mappings are not available.

When adding a new ingest:
1. Add download configuration to `download.yaml`
2. Create `src/<ingest_name>.py` with transform code
3. Create `src/<ingest_name>.yaml` with koza configuration
4. Add `<ingest_name>` to TRANSFORMS list in justfile
5. Create tests in `tests/test_<ingest_name>.py`
6. Update `src/versions.py` to declare the upstream source(s) and how to fetch their version

## Release Metadata

Every kozahub ingest emits an `output/release-metadata.yaml` describing the upstream sources, their versions, the artifacts produced, and the versions of build-time tools. This file is the contract monarch-ingest reads to assemble the merged knowledge graph's release receipt.

`src/versions.py` is the only per-ingest piece — it implements `get_source_versions()` returning a list of SourceVersion dicts. The `kozahub_metadata_schema` package provides reusable fetchers for the common patterns (HTTP Last-Modified, GitHub releases, URL-path regex, file-header parsing). The boilerplate (transform-content hashing, tool versions, build_version composition, yaml emission) is handled by `scripts/write_metadata.py`.

The `kozahub-metadata-schema` repo is expected as a sibling checkout (path-dep). Switch to a git or PyPI dep once published.

## Skills

- `.claude/skills/create-koza-ingest.md` - Create new koza ingests
- `.claude/skills/update-template.md` - Update to latest template version
