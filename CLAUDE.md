# pantherdb-orthologs-ingest

This is a Koza ingest repository for transforming PantherDB ortholog data into Biolink model format.

## Project Structure

- `download.yaml` - Configuration for downloading PantherDB ortholog data
- `src/` - Transform code and configuration
  - `transform.py` / `transform.yaml` - Main transform for ortholog associations
  - `panther_orthologs_utils.py` - Utility functions for parsing gene info and mapping
- `tests/` - Unit tests for transforms
- `output/` - Generated nodes and edges (gitignored)
- `data/` - Downloaded source data (gitignored)

## Key Commands

- `just run` - Full pipeline (download -> transform)
- `just download` - Download PantherDB ortholog data
- `just transform-all` - Run all transforms
- `just test` - Run tests

## Data Processing

This ingest processes the PantherDB AllOrthologs file, extracting orthology relationships between genes across multiple species. It uses an NCBI gene_info file to map species-specific identifiers to NCBIGene IDs when direct database mappings are not available.

When adding a new ingest:
1. Add download configuration to `download.yaml`
2. Create `src/<ingest_name>.py` with transform code
3. Create `src/<ingest_name>.yaml` with koza configuration
4. Add `<ingest_name>` to TRANSFORMS list in justfile
5. Create tests in `tests/test_<ingest_name>.py`

## Skills

- `.claude/skills/create-koza-ingest.md` - Create new koza ingests
- `.claude/skills/update-template.md` - Update to latest template version
