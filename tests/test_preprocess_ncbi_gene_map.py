"""
Tests for the NCBI gene map preprocessing SQL script.

These tests verify that the DuckDB-based preprocessing produces correct mappings
from composite keys (taxon_id|identifier) to NCBIGene IDs.
"""

import pytest
import duckdb
import gzip
import tempfile
import os
from pathlib import Path


@pytest.fixture
def test_gene_info_tsv(tmp_path):
    """Create a small test gene_info TSV with known values."""
    content = """#tax_id	GeneID	Symbol	LocusTag	Synonyms	dbXrefs	chromosome	map_location	description	type_of_gene	Symbol_from_nomenclature_authority	Full_name_from_nomenclature_authority	Nomenclature_status	Other_designations	Modification_date	Feature_type
227321	2868830	-	ANIA_08553	-	Ensembl:ANIA_08553	-	-	hypothetical protein	protein-coding	-	-	-	-	20240101	-
227321	2875778	AN0004	ANIA_00004	AN0004	Ensembl:ANIA_00004	-	-	some protein	protein-coding	-	-	-	-	20240101	-
9606	672	BRCA1	-	RNF53|BRCC1	MIM:113705|HGNC:HGNC:1100|Ensembl:ENSG00000012048	17	17q21.31	BRCA1 DNA repair	protein-coding	BRCA1	BRCA1 DNA repair associated	O	breast cancer type 1 susceptibility protein	20240101	-
10090	12189	Brca1	-	Brca1/Brca2-containing complex subunit 1	MGI:MGI:104537|Ensembl:ENSMUSG00000017146	11	11 B1.3	breast cancer 1	protein-coding	Brca1	breast cancer 1, early onset	O	breast cancer type 1 susceptibility protein homolog	20240101	-
"""
    gene_info_path = tmp_path / "test_gene_info.tsv"
    gene_info_path.write_text(content)
    return str(gene_info_path)


@pytest.fixture
def sql_template():
    """Load and parameterize the SQL template for testing."""
    sql_path = Path(__file__).parent.parent / "scripts" / "preprocess_ncbi_gene_map.sql"
    return sql_path.read_text()


def test_sql_produces_expected_mappings(test_gene_info_tsv, sql_template, tmp_path):
    """Test that the SQL produces correct composite key -> NCBIGene mappings."""
    # Modify SQL to use test file and output
    output_path = tmp_path / "test_output.tsv"

    # Adjust SQL for test file (uncompressed TSV)
    modified_sql = sql_template.replace(
        "'data/gene_info.gz'",
        f"'{test_gene_info_tsv}'"
    ).replace(
        "'data/ncbi_gene_map.tsv'",
        f"'{output_path}'"
    )

    # Execute the SQL
    con = duckdb.connect(":memory:")
    con.execute(modified_sql)

    # Read the output
    result = con.execute(f"SELECT * FROM read_csv('{output_path}', delim='\\t', header=true)").fetchall()
    mappings = {row[0]: row[1] for row in result}

    # Verify expected mappings exist
    # Aspergillus (227321) - ANIA_08553 should map to NCBIGene:2868830
    assert mappings.get("227321|ANIA_08553") == "NCBIGene:2868830"
    assert mappings.get("227321|AN0004") == "NCBIGene:2875778"
    assert mappings.get("227321|ANIA_00004") == "NCBIGene:2875778"

    # Human (9606) - BRCA1 mappings
    assert mappings.get("9606|BRCA1") == "NCBIGene:672"
    assert mappings.get("9606|RNF53") == "NCBIGene:672"
    assert mappings.get("9606|BRCC1") == "NCBIGene:672"
    assert mappings.get("9606|ENSG00000012048") == "NCBIGene:672"

    # Mouse (10090) - Brca1 mappings
    assert mappings.get("10090|Brca1") == "NCBIGene:12189"
    assert mappings.get("10090|ENSMUSG00000017146") == "NCBIGene:12189"


def test_sql_removes_ambiguous_keys(tmp_path):
    """Test that keys mapping to multiple genes are removed."""
    # Create test data where 'AmbiguousSymbol' maps to two different genes
    content = """#tax_id	GeneID	Symbol	LocusTag	Synonyms	dbXrefs	chromosome	map_location	description	type_of_gene	Symbol_from_nomenclature_authority	Full_name_from_nomenclature_authority	Nomenclature_status	Other_designations	Modification_date	Feature_type
9606	111	AmbiguousSymbol	-	-	-	-	-	gene 1	protein-coding	-	-	-	-	20240101	-
9606	222	AmbiguousSymbol	-	-	-	-	-	gene 2	protein-coding	-	-	-	-	20240101	-
9606	333	UniqueSymbol	-	-	-	-	-	gene 3	protein-coding	-	-	-	-	20240101	-
"""
    gene_info_path = tmp_path / "ambiguous_gene_info.tsv"
    gene_info_path.write_text(content)

    output_path = tmp_path / "ambiguous_output.tsv"

    sql_path = Path(__file__).parent.parent / "scripts" / "preprocess_ncbi_gene_map.sql"
    sql_template = sql_path.read_text()

    modified_sql = sql_template.replace(
        "'data/gene_info.gz'",
        f"'{gene_info_path}'"
    ).replace(
        "'data/ncbi_gene_map.tsv'",
        f"'{output_path}'"
    )

    con = duckdb.connect(":memory:")
    con.execute(modified_sql)

    result = con.execute(f"SELECT * FROM read_csv('{output_path}', delim='\\t', header=true)").fetchall()
    mappings = {row[0]: row[1] for row in result}

    # Ambiguous symbol should NOT be in the map
    assert "9606|AmbiguousSymbol" not in mappings

    # Unique symbol should still be present
    assert mappings.get("9606|UniqueSymbol") == "NCBIGene:333"


def test_sql_filters_irrelevant_taxons(tmp_path):
    """Test that only relevant taxons are included in the output."""
    # Create test data with an irrelevant taxon (99999)
    content = """#tax_id	GeneID	Symbol	LocusTag	Synonyms	dbXrefs	chromosome	map_location	description	type_of_gene	Symbol_from_nomenclature_authority	Full_name_from_nomenclature_authority	Nomenclature_status	Other_designations	Modification_date	Feature_type
99999	12345	IrrelevantGene	-	-	-	-	-	should be filtered	protein-coding	-	-	-	-	20240101	-
9606	54321	RelevantGene	-	-	-	-	-	should be included	protein-coding	-	-	-	-	20240101	-
"""
    gene_info_path = tmp_path / "taxon_filter_gene_info.tsv"
    gene_info_path.write_text(content)

    output_path = tmp_path / "taxon_filter_output.tsv"

    sql_path = Path(__file__).parent.parent / "scripts" / "preprocess_ncbi_gene_map.sql"
    sql_template = sql_path.read_text()

    modified_sql = sql_template.replace(
        "'data/gene_info.gz'",
        f"'{gene_info_path}'"
    ).replace(
        "'data/ncbi_gene_map.tsv'",
        f"'{output_path}'"
    )

    con = duckdb.connect(":memory:")
    con.execute(modified_sql)

    result = con.execute(f"SELECT * FROM read_csv('{output_path}', delim='\\t', header=true)").fetchall()
    mappings = {row[0]: row[1] for row in result}

    # Irrelevant taxon should NOT be in the map
    assert "99999|IrrelevantGene" not in mappings

    # Relevant taxon should be present
    assert mappings.get("9606|RelevantGene") == "NCBIGene:54321"
