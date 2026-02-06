"""
Unit tests for Panther Gene Orthology relationships ingest.

Tests the transform logic and parse_gene_info function using mock koza transforms.
"""

import gzip
from collections import Counter

import pytest
from biolink_model.datamodel.pydanticmodel_v2 import AgentTypeEnum, GeneToGeneHomologyAssociation, KnowledgeLevelEnum

from panther_orthologs_utils import (
    db_to_curie_map,
    panther_taxon_map,
    parse_gene_info,
)

# Columns needed for building test gene map (matches original preprocessing logic)
RELEVANT_NCBI_COLS = [
    "#tax_id",
    "GeneID",
    "Symbol",
    "LocusTag",
    "Synonyms",
    "dbXrefs",
    "Symbol_from_nomenclature_authority",
    "Full_name_from_nomenclature_authority",
    "Other_designations",
]

# Taxon IDs we care about (derived from panther_taxon_map values)
RELEVANT_NCBI_TAXONS = {v: "" for v in panther_taxon_map.values()}


def make_ncbi_taxon_gene_map(gene_info_file: str, relevant_columns: list, taxon_catalog: dict):
    """
    Build a nested dict mapping {taxon_id: {identifier: ncbi_gene_id}} from NCBI gene_info file.

    This function is used for testing only. In production, the preprocessing SQL script
    creates a flat TSV file that koza loads as a map.
    """
    # Ensure relevant columns has #tx_id as the first entry
    if relevant_columns[0] != "#tax_id":
        raise RuntimeError("- '#tax_id' must be first element present in relevant_columns arg... Exiting")

    # Many-->1 mapping dictionary
    taxa_gene_map = {tx_id: {} for tx_id in taxon_catalog}
    # Removes unreliable mapping keys from taxa_gene_map
    taxa_remove_map = {tx_id: Counter() for tx_id in taxon_catalog}

    with gzip.open(gene_info_file, "rt") as infile:
        # Read header line into memory to index relevant column fields
        hinfo = {hfield: i for i, hfield in enumerate(infile.readline().strip("\r").strip("\n").split("\t"))}

        # Now loop through each line, and create a map back to taxon / NCBIGene:xyz ...
        for line in infile:
            cols = line.strip("\r").strip("\n").split("\t")
            rel_data = [str(cols[hinfo[r]]) for r in relevant_columns]
            tx_id = rel_data[0]
            ncbi_gene_id = cols[hinfo["GeneID"]]

            # Only consume species we are interested in
            if tx_id not in taxon_catalog:
                continue

            # Find reliable mapping keys to this NCBI gene id
            # We take the set() of relevant mapping keys here... that if the same_id is reported on the same line
            # This removes the possibility of removing an id that is reported twice on the same line
            for map_key in set(rel_data[1:]):
                # Some columns like dbxref contain this character,
                # which separates common names from each other. So we loop through them
                # (minority, not majority have this in them)
                mk_cols = map_key.split("|")
                for key_to_ncbi in mk_cols:
                    # Deal with entries like MGI:MGI:95886, where we want to remove one of the MGI: prefix
                    key_split = key_to_ncbi.split(":")
                    if len(key_split) >= 2:
                        if key_split[0] == key_split[1]:
                            key_to_ncbi = "{}:{}".format(key_split[0], key_split[-1])

                    if key_to_ncbi not in taxa_gene_map[tx_id]:
                        taxa_gene_map[tx_id].update({key_to_ncbi: ncbi_gene_id})
                    else:
                        taxa_remove_map[tx_id][key_to_ncbi] += 1

    # Remove unreliable mapping keys to this NCBI gene id
    for tx_id in taxa_remove_map:
        for remove_key, rcount in taxa_remove_map[tx_id].items():
            del taxa_gene_map[tx_id][remove_key]

    # Return cleaned map back to a ncbi gene id that can be normalized later down road
    return taxa_gene_map


class MockKozaTransform:
    """Mock KozaTransform for testing that uses a fallback_map."""

    def __init__(self, fallback_map):
        self._fallback_map = fallback_map

    def lookup(self, key, value_column, map_name=None):
        """Lookup a value from the fallback map using composite key."""
        # Key format: "taxon_id|identifier"
        parts = key.split("|", 1)
        if len(parts) == 2:
            taxon_id, identifier = parts
            if taxon_id in self._fallback_map and identifier in self._fallback_map[taxon_id]:
                return f"NCBIGene:{self._fallback_map[taxon_id][identifier]}"
        # Return key if not found (matches koza behavior)
        return key


def run_transform(rows: list[dict], map_cache: dict) -> list:
    """Run the transform on a list of rows and return the results."""
    # Import parse_gene_info directly to test with mock koza transform
    from panther_orthologs_utils import db_to_curie_map, panther_taxon_map, parse_gene_info

    # Create a mock koza transform that uses the test map cache
    mock_koza = MockKozaTransform(map_cache)

    # Directly test the transform logic
    results = []
    for row in rows:
        # Replicate transform logic using the parse_gene_info function
        species_a, gene_a = parse_gene_info(row["Gene"], panther_taxon_map, db_to_curie_map, koza_transform=mock_koza)
        species_b, gene_b = parse_gene_info(
            row["Ortholog"], panther_taxon_map, db_to_curie_map, koza_transform=mock_koza
        )

        if (not species_a) or (not species_b):
            continue

        panther_ortholog_id = row["Panther Ortholog ID"]

        association = GeneToGeneHomologyAssociation(
            id="uuid:test",
            subject=gene_a,
            object=gene_b,
            predicate="biolink:orthologous_to",
            has_evidence=[f"PANTHER.FAMILY:{panther_ortholog_id}"],
            aggregator_knowledge_source=["infores:monarchinitiative"],
            primary_knowledge_source="infores:panther",
            knowledge_level=KnowledgeLevelEnum.knowledge_assertion,
            agent_type=AgentTypeEnum.not_provided,
        )
        results.append(association)
    return results


############################################################################
### Fixtures referencing code and basic parameters reused in tests below ###
@pytest.fixture()
def relevant_association_test_keys():
    """
    :return: list of keys that are relevant to the association test
    """
    return [
        "subject",
        "object",
        "predicate",
        "has_evidence",
        "aggregator_knowledge_source",
        "primary_knowledge_source",
        "knowledge_level",
        "agent_type",
    ]


@pytest.fixture
def map_cache():
    """Build a test gene map from the small test gene_info file."""
    return make_ncbi_taxon_gene_map(
        gene_info_file="./tests/test_ncbi_gene_info.txt.gz",
        relevant_columns=RELEVANT_NCBI_COLS,
        taxon_catalog=RELEVANT_NCBI_TAXONS,
    )


#############################################################################
### Fixture for panther rows/records to test for proper koza associations ###
@pytest.fixture
def panther_rows():
    data = [  # Human and rat ortholog row test
        (
            {
                "Gene": "HUMAN|HGNC=11477|UniProtKB=Q6GZX4",  # species1|DB=id1|protdb=pdbid1
                "Ortholog": "RAT|RGD=1564893|UniProtKB=Q6GZX2",  # species2|DB=id2|protdb=pdbid2
                "Type of ortholog": "LDO",  # [LDO, O, P, X ,LDX]  see: localtt
                "Common ancestor for the orthologs": "Euarchontoglires",  # unused
                "Panther Ortholog ID": "PTHR12434",
            },  # panther_id
            {
                "subject": "HGNC:11477",
                "object": "RGD:1564893",
                "predicate": "biolink:orthologous_to",
                "has_evidence": ["PANTHER.FAMILY:PTHR12434"],
                "aggregator_knowledge_source": ["infores:monarchinitiative"],
                "primary_knowledge_source": "infores:panther",
                "knowledge_level": KnowledgeLevelEnum.knowledge_assertion,
                "agent_type": AgentTypeEnum.not_provided,
            },
        ),
        # Mouse and Schizosaccharomyces pombe ortholog row test
        (
            {
                "Gene": "MOUSE|MGI=MGI=2147627|UniProtKB=Q91WQ3",
                "Ortholog": "SCHPO|PomBase=SPAC30C2.04|UniProtKB=Q9P6K7",
                "Type of ortholog": "LDO",
                "Common ancestor for the orthologs": "Opisthokonts",
                "Panther Ortholog ID": "PTHR11586",
            },
            {
                "subject": "MGI:2147627",
                "object": "PomBase:SPAC30C2.04",
                "predicate": "biolink:orthologous_to",
                "has_evidence": ["PANTHER.FAMILY:PTHR11586"],
                "aggregator_knowledge_source": ["infores:monarchinitiative"],
                "primary_knowledge_source": "infores:panther",
                "knowledge_level": KnowledgeLevelEnum.knowledge_assertion,
                "agent_type": AgentTypeEnum.not_provided,
            },
        ),
        # Xenopus tropicalis and ceravisea ("yeast") ortholog row test
        (
            {
                "Gene": "XENTR|Xenbase=XB-GENE-957143|UniProtKB=Q6P335",
                "Ortholog": "YEAST|SGD=S000004439|UniProtKB=P32366",
                "Type of ortholog": "O",
                "Common ancestor for the orthologs": "Opisthokonts",
                "Panther Ortholog ID": "PTHR11028",
            },
            {
                "subject": "Xenbase:XB-GENE-957143",
                "object": "SGD:S000004439",
                "predicate": "biolink:orthologous_to",
                "has_evidence": ["PANTHER.FAMILY:PTHR11028"],
                "aggregator_knowledge_source": ["infores:monarchinitiative"],
                "primary_knowledge_source": "infores:panther",
                "knowledge_level": KnowledgeLevelEnum.knowledge_assertion,
                "agent_type": AgentTypeEnum.not_provided,
            },
        ),
        # Zebrafish and fly ortholog row test
        (
            {
                "Gene": "DANRE|ZFIN=ZDB-GENE-050417-421|UniProtKB=Q567X8",
                "Ortholog": "DROME|FlyBase=FBgn0002773|UniProtKB=P18432",
                "Type of ortholog": "O",
                "Common ancestor for the orthologs": "Bilateria",
                "Panther Ortholog ID": "PTHR23049",
            },
            {
                "subject": "ZFIN:ZDB-GENE-050417-421",
                "object": "FB:FBgn0002773",
                "predicate": "biolink:orthologous_to",
                "has_evidence": ["PANTHER.FAMILY:PTHR23049"],
                "aggregator_knowledge_source": ["infores:monarchinitiative"],
                "primary_knowledge_source": "infores:panther",
                "knowledge_level": KnowledgeLevelEnum.knowledge_assertion,
                "agent_type": AgentTypeEnum.not_provided,
            },
        ),
        # C. elegans and Dictyostelium discoideum ortholog row test
        (
            {
                "Gene": "CAEEL|WormBase=WBGene00009059|UniProtKB=Q19739",
                "Ortholog": "DICDI|dictyBase=DDB_G0269178|UniProtKB=Q9GPS0",
                "Type of ortholog": "O",
                "Common ancestor for the orthologs": "Unikonts",
                "Panther Ortholog ID": "PTHR24072",
            },
            {
                "subject": "WB:WBGene00009059",
                "object": "dictyBase:DDB_G0269178",
                "predicate": "biolink:orthologous_to",
                "has_evidence": ["PANTHER.FAMILY:PTHR24072"],
                "aggregator_knowledge_source": ["infores:monarchinitiative"],
                "primary_knowledge_source": "infores:panther",
                "knowledge_level": KnowledgeLevelEnum.knowledge_assertion,
                "agent_type": AgentTypeEnum.not_provided,
            },
        ),
        # Human Ensembl gene and mouse ortholog row test
        (
            {
                "Gene": "HUMAN|Ensembl=ENSG00000275949.5|UniProtKB=A0A0G2JMH3",
                "Ortholog": "MOUSE|MGI=MGI=99431|UniProtKB=P84078",
                "Type of ortholog": "O",
                "Common ancestor for the orthologs": "Euarchontoglires",
                "Panther Ortholog ID": "PTHR11711",
            },
            {
                "subject": "ENSEMBL:ENSG00000275949",
                "object": "MGI:99431",
                "predicate": "biolink:orthologous_to",
                "has_evidence": ["PANTHER.FAMILY:PTHR11711"],
                "aggregator_knowledge_source": ["infores:monarchinitiative"],
                "primary_knowledge_source": "infores:panther",
                "knowledge_level": KnowledgeLevelEnum.knowledge_assertion,
                "agent_type": AgentTypeEnum.not_provided,
            },
        ),
    ]

    return data


#################################################
### Fixtures to test parse_gene_info function ###
@pytest.fixture
def panther_species_genes():
    species_genes = [
        ["HUMAN|Ensembl=ENSG00000275949.5|UniProtKB=A0A0G2JMH3", {"ENSEMBL:ENSG00000275949": ""}],
        ["HUMAN|HGNC=11477|UniProtKB=Q6GZX4", {"HGNC:11477": ""}],
        ["MOUSE|MGI=MGI=99431|UniProtKB=P84078", {"MGI:99431": ""}],
        ["CANLF|Ensembl=ENSCAFG00845004769.1", {"ENSEMBL:ENSCAFG00845004769": ""}],
        ["BOVIN|Ensembl=ENSBTAG00000048390.1|UniProtKB=A0A3Q1NMJ3", {"ENSEMBL:ENSBTAG00000048390": ""}],
        ["PIG|Ensembl=ENSSSCG00000033574.3|UniProtKB=A0A8W4FPJ3", {"ENSEMBL:ENSSSCG00000033574": ""}],
        ["RAT|Ensembl=ENSRNOG00000066524.1|UniProtKB=A0A8I5ZRQ5", {"ENSEMBL:ENSRNOG00000066524": ""}],
        ["CHICK|Ensembl=ENSGALG00000014680|UniProtKB=F1NB96", {"ENSEMBL:ENSGALG00000014680": ""}],
        ["XENTR|Ensembl=ENSXETG00000030579|UniProtKB=A0A6I8PUG3", {"ENSEMBL:ENSXETG00000030579": ""}],
        ["DANRE|ZFIN=ZDB-GENE-080205-1|UniProtKB=A8WFS6", {"ZFIN:ZDB-GENE-080205-1": ""}],
        ["DROME|FlyBase=FBgn0010348|UniProtKB=P61209", {"FB:FBgn0010348": ""}],
        ["CAEEL|WormBase=WBGene00000446|UniProtKB=P34663", {"WB:WBGene00000446": ""}],
        ["DICDI|dictyBase=DDB_G0274381|UniProtKB=P54642", {"dictyBase:DDB_G0274381": ""}],
        # Aspergillus will be our test cases to ensure mapping back to NCBIGene is done properly
        # and scenario where we use UniProtKB identifier instead
        ["EMENI|Gene_ORFName=AN0062|UniProtKB=Q5BHB8", {"NCBIGene:ANIA_00062": "", "UniProtKB:Q5BHB8": ""}],
        ["EMENI|EnsemblGenome=ANIA_08553|UniProtKB=Q5AT27", {"NCBIGene:2868830": "", "UniProtKB:Q5AT27": ""}],
        ["SCHPO|PomBase=SPAC13G7.02c|UniProtKB=Q10265", {"PomBase:SPAC13G7.02c": ""}],
        ["YEAST|SGD=S000003465|UniProtKB=P17442", {"SGD:S000003465": ""}],
    ]

    return species_genes


@pytest.fixture
def exclude_species_genes():
    exclude_species_genes = ["FAKE_SPECIES_13|Ensembl=ENSFAKE00000000001|UniProtKB=FAKE1234"]
    return exclude_species_genes


###########################################################
### Perform our actual tests here on the fixtures above ###
def test_panther_rows(panther_rows, relevant_association_test_keys, map_cache):
    # Upack our test rows and expected info
    rows_to_test, expected_res = [v[0] for v in panther_rows], [v[1] for v in panther_rows]

    # Run the transform using the new KozaRunner pattern
    koza_associations = run_transform(rows_to_test, map_cache)

    # Now check our koza generated associations against expected results
    for koza_association, expected_info in zip(koza_associations, expected_res):
        # Ensure association type is correct (we are dealing with GeneToGeneHomologyAssociation s here)
        assert isinstance(koza_association, GeneToGeneHomologyAssociation)

        # Test that koza processing of pantherdb row produces expected results
        for key in relevant_association_test_keys:
            assert key in expected_info
            assert getattr(koza_association, key) == expected_info[key]


def test_species_parse_gene(panther_species_genes, map_cache):
    # Create mock koza transform for testing
    mock_koza = MockKozaTransform(map_cache)

    for gene_info, expected in panther_species_genes:
        species, gene_id = parse_gene_info(gene_info, panther_taxon_map, db_to_curie_map, koza_transform=mock_koza)

        # Assert that the parsed species and gene ID match the expected values
        assert species in panther_taxon_map
        assert gene_id in expected  # Allows for multiple mapping values to be present


def test_exclude_species_parse_gene(exclude_species_genes, map_cache):
    # Create mock koza transform for testing
    mock_koza = MockKozaTransform(map_cache)

    for gene_info in exclude_species_genes:
        species, gene_id = parse_gene_info(gene_info, panther_taxon_map, db_to_curie_map, koza_transform=mock_koza)

        # Assert that the species and gene ID are empty for excluded species
        assert not species
        assert not gene_id
