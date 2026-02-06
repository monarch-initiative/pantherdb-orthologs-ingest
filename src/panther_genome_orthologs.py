"""
Ingest of Reference Genome Orthologs from Panther
"""

# Imports
import sys
from pathlib import Path

# Add the src directory to sys.path for relative imports
_src_dir = Path(__file__).parent
if str(_src_dir) not in sys.path:
    sys.path.insert(0, str(_src_dir))

import uuid  # noqa: E402

import koza  # noqa: E402
from biolink_model.datamodel.pydanticmodel_v2 import (  # noqa: E402
    AgentTypeEnum,
    GeneToGeneHomologyAssociation,
    KnowledgeLevelEnum,
)

# Custom pantherdb specific function, and constants respectively
from panther_orthologs_utils import (  # noqa: E402
    db_to_curie_map,
    panther_taxon_map,
    parse_gene_info,
)


@koza.transform_record()
def transform(koza_transform, row: dict) -> list[GeneToGeneHomologyAssociation]:
    """Transform a single row from the Panther orthologs file."""
    # Parse the gene information for both species and format gene id to curie:gene_id
    # (Gene and Ortholog columns are formatted the same, but for different species/gene info)
    # Uses koza_transform.lookup() for NCBI gene map lookups
    species_a, gene_a = parse_gene_info(row["Gene"], panther_taxon_map, db_to_curie_map, koza_transform=koza_transform)
    species_b, gene_b = parse_gene_info(
        row["Ortholog"], panther_taxon_map, db_to_curie_map, koza_transform=koza_transform
    )

    # Only consume species we are interested in (i.e., those that are in our ncbitaxon_catalog)
    if (not species_a) or (not species_b):
        return []

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
        agent_type=AgentTypeEnum.not_provided,
    )

    return [association]
