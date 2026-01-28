"""
Ingest of Reference Genome Orthologs from Panther
"""

# Imports
import os
import uuid
import koza
from biolink_model.datamodel.pydanticmodel_v2 import GeneToGeneHomologyAssociation, KnowledgeLevelEnum, AgentTypeEnum

# Custom pantherdb specific function, and constants respectively
from panther_orthologs_utils import (make_ncbi_taxon_gene_map, parse_gene_info)
from panther_orthologs_utils import (panther_taxon_map, relevant_ncbi_cols, relevant_ncbi_taxons, db_to_curie_map)


# Lazy-load the NCBI gene map to avoid loading during test imports
_tx_gmap = None


def get_tx_gmap():
    """Get the taxon gene map, loading it if not already loaded."""
    global _tx_gmap
    if _tx_gmap is None:
        ncbi_map_file = "./data/gene_info.gz"
        _tx_gmap = make_ncbi_taxon_gene_map(
            gene_info_file=ncbi_map_file,
            relevant_columns=relevant_ncbi_cols,
            taxon_catalog=relevant_ncbi_taxons
        )
    return _tx_gmap


# For testing: allow injection of a custom map
def set_tx_gmap(gmap):
    """Set the taxon gene map (used for testing)."""
    global _tx_gmap
    _tx_gmap = gmap


@koza.transform_record()
def transform(koza_transform, row: dict) -> GeneToGeneHomologyAssociation | None:
    """Transform a single row from the Panther orthologs file."""
    tx_gmap = get_tx_gmap()

    # Parse the gene information for both species and format gene id to curie:gene_id
    # (Gene and Ortholog columns are formatted the same, but for different species/gene info)
    species_a, gene_a = parse_gene_info(row["Gene"], panther_taxon_map, db_to_curie_map, tx_gmap)
    species_b, gene_b = parse_gene_info(row["Ortholog"], panther_taxon_map, db_to_curie_map, tx_gmap)

    # Only consume species we are interested in (i.e., those that are in our ncbitaxon_catalog)
    if (not species_a) or (not species_b):
        return None

    # Our ortholog identifier (panther protein family name), and predicate
    panther_ortholog_id = row["Panther Ortholog ID"]
    predicate = "biolink:orthologous_to"

    # Generate our association object (uuid4 for reliable uniqueness across environments)
    association = GeneToGeneHomologyAssociation(
        id="uuid:{}".format(str(uuid.uuid4())),
        subject=gene_a,
        object=gene_b,
        predicate=predicate,
        has_evidence=["PANTHER.FAMILY:{}".format(panther_ortholog_id)],
        aggregator_knowledge_source=["infores:monarchinitiative"],
        primary_knowledge_source="infores:panther",
        knowledge_level=KnowledgeLevelEnum.knowledge_assertion,
        agent_type=AgentTypeEnum.not_provided
    )

    return association