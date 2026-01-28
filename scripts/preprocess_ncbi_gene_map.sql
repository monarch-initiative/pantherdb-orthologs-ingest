-- Preprocessing SQL for NCBI gene_info to create a flat gene map
-- Input: data/gene_info.gz (tab-separated, with header)
-- Output: data/ncbi_gene_map.tsv (composite_key -> ncbi_gene_id)
--
-- The map allows lookups of NCBIGene IDs using a composite key of taxon|identifier
-- where identifier can be Symbol, LocusTag, or a parsed dbXrefs value.

-- Relevant taxon IDs for the pantherdb orthologs ingest
-- These correspond to species in the panther_taxon_map
CREATE TEMP TABLE relevant_taxons AS
SELECT * FROM (VALUES
    (9606),   -- HUMAN
    (10090),  -- MOUSE
    (9615),   -- CANLF (dog)
    (9913),   -- BOVIN (cow)
    (9823),   -- PIG
    (10116),  -- RAT
    (9031),   -- CHICK
    (8364),   -- XENTR (Xenopus tropicalis)
    (7955),   -- DANRE (zebrafish)
    (7227),   -- DROME (fly)
    (6239),   -- CAEEL (C. elegans)
    (44689),  -- DICDI (Dictyostelium)
    (227321), -- EMENI (Aspergillus)
    (4896),   -- SCHPO (S. pombe)
    (4932)    -- YEAST (S. cerevisiae)
) AS t(tax_id);

-- Read gene_info with header (DuckDB will auto-detect the 16 columns)
-- We need all columns that the original Python function uses for mapping:
-- Symbol, LocusTag, Synonyms, dbXrefs, Symbol_from_nomenclature_authority,
-- Full_name_from_nomenclature_authority, Other_designations
CREATE TEMP TABLE gene_info_raw AS
SELECT
    CAST("#tax_id" AS INTEGER) AS tax_id,
    "GeneID" AS gene_id,
    "Symbol" AS symbol,
    "LocusTag" AS locus_tag,
    "Synonyms" AS synonyms,
    "dbXrefs" AS db_xrefs,
    "Symbol_from_nomenclature_authority" AS symbol_authority,
    "Full_name_from_nomenclature_authority" AS full_name_authority,
    "Other_designations" AS other_designations
FROM read_csv(
    'data/gene_info.gz',
    delim = '\t',
    header = true,
    auto_detect = true
);

-- Filter to relevant taxons only
CREATE TEMP TABLE gene_info AS
SELECT g.*
FROM gene_info_raw g
INNER JOIN relevant_taxons t ON g.tax_id = t.tax_id;

-- Create mappings from Symbol (if not '-')
CREATE TEMP TABLE symbol_mappings AS
SELECT
    tax_id || '|' || symbol AS composite_key,
    'NCBIGene:' || gene_id AS ncbi_gene_id
FROM gene_info
WHERE symbol != '-';

-- Create mappings from LocusTag (if not '-')
CREATE TEMP TABLE locus_tag_mappings AS
SELECT
    tax_id || '|' || locus_tag AS composite_key,
    'NCBIGene:' || gene_id AS ncbi_gene_id
FROM gene_info
WHERE locus_tag != '-';

-- Parse dbXrefs to extract additional identifiers
-- dbXrefs format: "DB1:ID1|DB2:ID2|..." or "-"
-- We want to create mappings for each part
CREATE TEMP TABLE dbxref_exploded AS
SELECT
    tax_id,
    gene_id,
    UNNEST(string_split(db_xrefs, '|')) AS db_xref_part
FROM gene_info
WHERE db_xrefs != '-';

-- Create mappings from parsed dbXrefs
-- We use the part after the colon as the identifier
CREATE TEMP TABLE dbxref_mappings AS
SELECT
    tax_id || '|' ||
    CASE
        -- Handle special case for Ensembl versioned IDs (ANIA_08553 instead of Ensembl:ANIA_08553)
        WHEN db_xref_part LIKE 'Ensembl%' THEN
            regexp_extract(db_xref_part, ':(.+)', 1)
        -- For Gene_ORFName like AN0062
        WHEN db_xref_part LIKE 'HPRD:%' THEN
            regexp_extract(db_xref_part, ':(.+)', 1)
        -- Default: use the full value after the colon
        ELSE
            regexp_extract(db_xref_part, ':(.+)', 1)
    END AS composite_key,
    'NCBIGene:' || gene_id AS ncbi_gene_id
FROM dbxref_exploded
WHERE db_xref_part LIKE '%:%'
  AND regexp_extract(db_xref_part, ':(.+)', 1) IS NOT NULL
  AND regexp_extract(db_xref_part, ':(.+)', 1) != '';

-- Create mappings from Synonyms (pipe-separated, if not '-')
CREATE TEMP TABLE synonyms_exploded AS
SELECT
    tax_id,
    gene_id,
    UNNEST(string_split(synonyms, '|')) AS synonym
FROM gene_info
WHERE synonyms != '-';

CREATE TEMP TABLE synonym_mappings AS
SELECT
    tax_id || '|' || synonym AS composite_key,
    'NCBIGene:' || gene_id AS ncbi_gene_id
FROM synonyms_exploded
WHERE synonym IS NOT NULL AND synonym != '';

-- Create mappings from Symbol_from_nomenclature_authority (if not '-')
CREATE TEMP TABLE symbol_authority_mappings AS
SELECT
    tax_id || '|' || symbol_authority AS composite_key,
    'NCBIGene:' || gene_id AS ncbi_gene_id
FROM gene_info
WHERE symbol_authority != '-';

-- Create mappings from Full_name_from_nomenclature_authority (if not '-')
CREATE TEMP TABLE full_name_authority_mappings AS
SELECT
    tax_id || '|' || full_name_authority AS composite_key,
    'NCBIGene:' || gene_id AS ncbi_gene_id
FROM gene_info
WHERE full_name_authority != '-';

-- Create mappings from Other_designations (pipe-separated, if not '-')
CREATE TEMP TABLE other_designations_exploded AS
SELECT
    tax_id,
    gene_id,
    UNNEST(string_split(other_designations, '|')) AS designation
FROM gene_info
WHERE other_designations != '-';

CREATE TEMP TABLE other_designation_mappings AS
SELECT
    tax_id || '|' || designation AS composite_key,
    'NCBIGene:' || gene_id AS ncbi_gene_id
FROM other_designations_exploded
WHERE designation IS NOT NULL AND designation != '';

-- Combine all mappings
-- All sources are combined: Symbol, LocusTag, dbXrefs, Synonyms,
-- Symbol_from_nomenclature_authority, Full_name_from_nomenclature_authority, Other_designations
CREATE TEMP TABLE all_mappings AS
SELECT DISTINCT composite_key, ncbi_gene_id
FROM (
    SELECT composite_key, ncbi_gene_id FROM symbol_mappings
    UNION ALL
    SELECT composite_key, ncbi_gene_id FROM locus_tag_mappings
    UNION ALL
    SELECT composite_key, ncbi_gene_id FROM dbxref_mappings
    UNION ALL
    SELECT composite_key, ncbi_gene_id FROM synonym_mappings
    UNION ALL
    SELECT composite_key, ncbi_gene_id FROM symbol_authority_mappings
    UNION ALL
    SELECT composite_key, ncbi_gene_id FROM full_name_authority_mappings
    UNION ALL
    SELECT composite_key, ncbi_gene_id FROM other_designation_mappings
) combined
WHERE composite_key IS NOT NULL
  AND composite_key NOT LIKE '%|'
  AND composite_key NOT LIKE '|%';

-- Find ambiguous keys (same composite_key maps to different gene IDs)
-- These should be removed for reliability, matching the original Python logic
CREATE TEMP TABLE ambiguous_keys AS
SELECT composite_key
FROM all_mappings
GROUP BY composite_key
HAVING COUNT(DISTINCT ncbi_gene_id) > 1;

-- Final map: remove ambiguous keys
CREATE TABLE ncbi_gene_map AS
SELECT m.composite_key, m.ncbi_gene_id
FROM all_mappings m
LEFT JOIN ambiguous_keys a ON m.composite_key = a.composite_key
WHERE a.composite_key IS NULL
ORDER BY m.composite_key;

-- Export to TSV
COPY ncbi_gene_map TO 'data/ncbi_gene_map.tsv' (DELIMITER '\t', HEADER);
