# pantherdb-orthologs-ingest justfile

# Package directory
PKG := "src"

# Explicitly enumerate transforms
TRANSFORMS := "panther_genome_orthologs"

# List all commands
_default:
    @just --list

# ============== Project Management ==============

# Install dependencies
[group('project management')]
install:
    uv sync --group dev

# ============== Ingest Pipeline ==============

# Full pipeline: download -> preprocess -> transform -> postprocess
[group('ingest')]
run: download preprocess transform-all postprocess
    @echo "Done!"

# Download source data using kghub-downloader
[group('ingest')]
download: install
    uv run downloader download.yaml

# Preprocess: create koza-compatible map files from gene_info.gz
[group('ingest')]
preprocess:
    uv run python scripts/preprocess_ncbi_gene_map.py

# Run all transforms
[group('ingest')]
transform-all:
    #!/usr/bin/env bash
    set -euo pipefail
    for t in {{TRANSFORMS}}; do
        if [ -n "$t" ]; then
            echo "Transforming $t..."
            uv run koza transform {{PKG}}/$t.yaml
        fi
    done

# Run specific transform
[group('ingest')]
transform NAME:
    uv run koza transform {{PKG}}/{{NAME}}.yaml

# Postprocess (no-op for this ingest)
[group('ingest')]
postprocess:
    @echo "No postprocessing required"

# ============== Development ==============

# Run tests
[group('development')]
test:
    uv run pytest

# Run tests with coverage
[group('development')]
test-cov:
    uv run pytest --cov=. --cov-report=term-missing

# Lint code
[group('development')]
lint:
    uv run ruff check .

# Format code
[group('development')]
format:
    uv run ruff format .

# Clean output directory
[group('development')]
clean:
    rm -rf output/
