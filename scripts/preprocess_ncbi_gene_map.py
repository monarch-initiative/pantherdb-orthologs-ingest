#!/usr/bin/env python
"""
Preprocess NCBI gene_info.gz to create a flat gene map for koza lookups.

This script runs the SQL preprocessing using DuckDB to create a TSV file
that maps composite keys (taxon|identifier) to NCBIGene IDs.

Input: data/gene_info.gz
Output: data/ncbi_gene_map.tsv
"""

from pathlib import Path

import duckdb


def main():
    """Run the preprocessing SQL script."""
    script_dir = Path(__file__).parent
    project_root = script_dir.parent
    sql_file = script_dir / "preprocess_ncbi_gene_map.sql"

    # Change to project root so relative paths in SQL work correctly
    import os
    original_cwd = os.getcwd()
    os.chdir(project_root)

    try:
        # Read and execute the SQL script
        sql_content = sql_file.read_text()

        print("Preprocessing NCBI gene_info.gz to create gene map...")
        print(f"Input: data/gene_info.gz")
        print(f"Output: data/ncbi_gene_map.tsv")

        con = duckdb.connect(":memory:")
        con.execute(sql_content)

        # Get stats from the result table
        result = con.execute("SELECT COUNT(*) FROM ncbi_gene_map").fetchone()
        print(f"Created {result[0]:,} mappings")

        con.close()
        print("Done!")

    finally:
        os.chdir(original_cwd)


if __name__ == "__main__":
    main()
