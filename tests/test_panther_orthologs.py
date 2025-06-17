"""
An example test file for the transform script.

It uses pytest fixtures to define the input data and the mock koza transform.
The test_example function then tests the output of the transform script.

See the Koza documentation for more information on testing transforms:
https://koza.monarchinitiative.org/Usage/testing/

Unit tests for Panther Gene Orthology relationships ingest
"""


import pytest
from biolink_model.datamodel.pydanticmodel_v2 import GeneToGeneHomologyAssociation, KnowledgeLevelEnum, AgentTypeEnum
from koza.utils.testing_utils import mock_koza
from src.pantherdb_orthologs_ingest.panther_ortholog_utils import (make_ncbi_taxon_gene_map, 
                                                                   parse_gene_info)

from src.pantherdb_orthologs_ingest.panther_ortholog_utils import (panther_taxon_map, 
                                                                   relevant_ncbi_cols, 
                                                                   relevant_ncbi_taxons, 
                                                                   db_to_curie_map)


############################################################################
### Fixtures referencing code and basic parameters reused in tests below ###
@pytest.fixture
def source_name():
    """
    :return: string name of ingest source found within transform yaml
    """
    return "pantherdb_orthologs"


@pytest.fixture
def script():
    """
    :return: string path to Panther Gene Orthology relationships ingest script
    """
    return "./src/pantherdb_orthologs_ingest/transform.py"


@pytest.fixture()
def relevant_association_test_keys():
    """
    :return: list of keys that are relevant to the association test
    """
    return ["subject",
            "object",
            "predicate",
            "has_evidence",
            "aggregator_knowledge_source",
            "primary_knowledge_source",
            "knowledge_level",
            "agent_type"]


###########################################################################################
### Records to test for proper associations created each species (two species per test) ###
@pytest.fixture
def panther_row_v1(mock_koza, source_name, script):
    """
    Tests a row with a human gene and a rat ortholog
    """
    row = {"Gene": "HUMAN|HGNC=11477|UniProtKB=Q6GZX4",  # species1|DB=id1|protdb=pdbid1
           "Ortholog": "RAT|RGD=1564893|UniProtKB=Q6GZX2",  # species2|DB=id2|protdb=pdbid2
           "Type of ortholog": "LDO",  # [LDO, O, P, X ,LDX]  see: localtt
           "Common ancestor for the orthologs": "Euarchontoglires",  # unused
           "Panther Ortholog ID": "PTHR12434"}  # panther_id
    
    expected = {"subject": "HGNC:11477",
                "object": "RGD:1564893",
                "predicate": "biolink:orthologous_to",
                "has_evidence": ["PANTHER.FAMILY:PTHR12434"],
                "aggregator_knowledge_source": ["infores:monarchinitiative"],
                "primary_knowledge_source": "infores:panther",
                "knowledge_level": KnowledgeLevelEnum.knowledge_assertion,
                "agent_type": AgentTypeEnum.not_provided}
           
    return mock_koza(name=source_name, data=row, transform_code=script), expected


@pytest.fixture
def panther_row_v2(mock_koza, source_name, script):
    """
    Tests a row with a mouse gene and a Schizosaccharomyces pombe ortholog
    """
    row = {"Gene":"MOUSE|MGI=MGI=2147627|UniProtKB=Q91WQ3",
           "Ortholog": "SCHPO|PomBase=SPAC30C2.04|UniProtKB=Q9P6K7",
           "Type of ortholog": "LDO", 
           "Common ancestor for the orthologs": "Opisthokonts",
           "Panther Ortholog ID": "PTHR11586"}
    
    expected = {"subject": "MGI:2147627",
                "object": "PomBase:SPAC30C2.04",
                "predicate": "biolink:orthologous_to",
                "has_evidence": ["PANTHER.FAMILY:PTHR11586"],
                "aggregator_knowledge_source": ["infores:monarchinitiative"],
                "primary_knowledge_source": "infores:panther",
                "knowledge_level": KnowledgeLevelEnum.knowledge_assertion,
                "agent_type": AgentTypeEnum.not_provided}
           
    return mock_koza(name=source_name, data=row, transform_code=script), expected


@pytest.fixture
def panther_row_v3(mock_koza, source_name, script):
    """
    Tests a row with a Xenopus tropicalis gene and a yeast ortholog
    """
    row = {"Gene": "XENTR|Xenbase=XB-GENE-957143|UniProtKB=Q6P335",
           "Ortholog": "YEAST|SGD=S000004439|UniProtKB=P32366",
           "Type of ortholog": "O", 
           "Common ancestor for the orthologs": "Opisthokonts",
           "Panther Ortholog ID": "PTHR11028"}
    
    expected = {"subject": "Xenbase:XB-GENE-957143",
                "object": "SGD:S000004439",
                "predicate": "biolink:orthologous_to",
                "has_evidence": ["PANTHER.FAMILY:PTHR11028"],
                "aggregator_knowledge_source": ["infores:monarchinitiative"],
                "primary_knowledge_source": "infores:panther",
                "knowledge_level": KnowledgeLevelEnum.knowledge_assertion,
                "agent_type": AgentTypeEnum.not_provided}
    
    return mock_koza(name=source_name, data=row, transform_code=script), expected


@pytest.fixture
def panther_row_v4(mock_koza, source_name, script):
    """
    Tests a row with a zebrafish gene and a fly ortholog
    """
    row = {"Gene": "DANRE|ZFIN=ZDB-GENE-050417-421|UniProtKB=Q567X8",
           "Ortholog": "DROME|FlyBase=FBgn0002773|UniProtKB=P18432",
           "Type of ortholog": "O", 
           "Common ancestor for the orthologs": "Bilateria",
           "Panther Ortholog ID": "PTHR23049"}
    
    expected = {"subject": "ZFIN:ZDB-GENE-050417-421",
                "object": "FB:FBgn0002773",
                "predicate": "biolink:orthologous_to",
                "has_evidence": ["PANTHER.FAMILY:PTHR23049"],
                "aggregator_knowledge_source": ["infores:monarchinitiative"],
                "primary_knowledge_source": "infores:panther",
                "knowledge_level": KnowledgeLevelEnum.knowledge_assertion,
                "agent_type": AgentTypeEnum.not_provided}
    
    return mock_koza(name=source_name, data=row, transform_code=script), expected


@pytest.fixture
def panther_row_v5(mock_koza, source_name, script):
    """
    Tests c. elegans gene and a Dictyostelium discoideum ortholog
    """
    row = {"Gene": "CAEEL|WormBase=WBGene00009059|UniProtKB=Q19739",
           "Ortholog": "DICDI|dictyBase=DDB_G0269178|UniProtKB=Q9GPS0",
           "Type of ortholog": "O", 
           "Common ancestor for the orthologs": "Unikonts",
           "Panther Ortholog ID": "PTHR24072"}
    
    expected = {"subject": "WB:WBGene00009059",
                "object": "dictyBase:DDB_G0269178",
                "predicate": "biolink:orthologous_to",
                "has_evidence": ["PANTHER.FAMILY:PTHR24072"],
                "aggregator_knowledge_source": ["infores:monarchinitiative"],
                "primary_knowledge_source": "infores:panther",
                "knowledge_level": KnowledgeLevelEnum.knowledge_assertion,
                "agent_type": AgentTypeEnum.not_provided}
    
    return mock_koza(name=source_name, data=row, transform_code=script), expected


@pytest.fixture
def panther_row_v6(mock_koza, source_name, script):
    """
    Tests Ensembl identifier and the removal of transcript version numbers from the end of the gene identifier.
    """

    row = {"Gene": "HUMAN|Ensembl=ENSG00000275949.5|UniProtKB=A0A0G2JMH3",
           "Ortholog": "MOUSE|MGI=MGI=99431|UniProtKB=P84078",
           "Type of ortholog": "O", 
           "Common ancestor for the orthologs": "Euarchontoglires",
           "Panther Ortholog ID": "PTHR11711"}
    
    expected = {"subject": "ENSEMBL:ENSG00000275949",
                "object": "MGI:99431",
                "predicate": "biolink:orthologous_to",
                "has_evidence": ["PANTHER.FAMILY:PTHR11711"],
                "aggregator_knowledge_source": ["infores:monarchinitiative"],
                "primary_knowledge_source": "infores:panther",
                "knowledge_level": KnowledgeLevelEnum.knowledge_assertion,
                "agent_type": AgentTypeEnum.not_provided}
    
    return mock_koza(name=source_name, data=row, transform_code=script), expected


@pytest.fixture
def panther_row_fixtures():
    """
    Fixture to hold all record fixtures for testing.
    """
    return [panther_row_v1,
            panther_row_v2,
            panther_row_v3,
            panther_row_v4,
            panther_row_v5,
            panther_row_v6]


#################################################
### Fixtures to test parse_gene_info function ###
@pytest.fixture
def panther_species_genes():
    species_genes = [["HUMAN|Ensembl=ENSG00000275949.5|UniProtKB=A0A0G2JMH3", {"ENSEMBL:ENSG00000275949":''}],
                     ["HUMAN|HGNC=11477|UniProtKB=Q6GZX4", {"HGNC:11477":''}],
                     ["MOUSE|MGI=MGI=99431|UniProtKB=P84078", {"MGI:99431":''}],
                     ["CANLF|Ensembl=ENSCAFG00845004769.1", {"ENSEMBL:ENSCAFG00845004769":''}],
                     ["BOVIN|Ensembl=ENSBTAG00000048390.1|UniProtKB=A0A3Q1NMJ3", {"ENSEMBL:ENSBTAG00000048390":''}],
                     ["PIG|Ensembl=ENSSSCG00000033574.3|UniProtKB=A0A8W4FPJ3", {"ENSEMBL:ENSSSCG00000033574":''}],
                     ["RAT|Ensembl=ENSRNOG00000066524.1|UniProtKB=A0A8I5ZRQ5", {"ENSEMBL:ENSRNOG00000066524":''}],
                     ["CHICK|Ensembl=ENSGALG00000014680|UniProtKB=F1NB96", {"ENSEMBL:ENSGALG00000014680":''}],
                     ["XENTR|Ensembl=ENSXETG00000030579|UniProtKB=A0A6I8PUG3", {"ENSEMBL:ENSXETG00000030579":''}],
                     ["DANRE|ZFIN=ZDB-GENE-080205-1|UniProtKB=A8WFS6", {"ZFIN:ZDB-GENE-080205-1":''}],
                     ["DROME|FlyBase=FBgn0010348|UniProtKB=P61209", {"FB:FBgn0010348":''}],
                     ["CAEEL|WormBase=WBGene00000446|UniProtKB=P34663", {"WormBase:WBGene00000446":''}],
                     ["DICDI|dictyBase=DDB_G0274381|UniProtKB=P54642", {"dictyBase:DDB_G0274381":''}],
                     ["EMENI|Gene_ORFName=AN0062|UniProtKB=Q5BHB8", ""], ##??? Need to figure this one out still...
                     ["SCHPO|PomBase=SPAC13G7.02c|UniProtKB=Q10265", {"PomBase:SPAC13G7.02c":''}],
                     ["YEAST|SGD=S000003465|UniProtKB=P17442", {"SGD:S000003465":''}]]

    return species_genes


@pytest.fixture
def exclude_species_genes():
    exclude_species_genes = ["FAKE_SPECIES_13|Ensembl=ENSFAKE00000000001|UniProtKB=FAKE1234"]
    return exclude_species_genes


###########################################################
### Perform our actual tests here on the fixtures above ###
def test_all_records(mock_koza, source_name, script, panther_row_fixtures, rel_keys=relevant_association_test_keys):

    # Loop through all of our record testing fixtures, and perform our koza tests
    for rec_fixture in panther_record_fixtures:

        # Grab the mock Koza application and expected info
        koza_association, expected_info = rec_fixture(mock_koza, source_name, script)

        # Ensure association type is correct (we are dealing with GeneToGeneHomologyAssociation s here)
        assert len(koza_assocation) == 1
        assert isinstance(koza_association[0], GeneToGeneHomologyAssociation)

        # Test that mock koza processing of pantherdb row produces expected results
        for key in rel_keys:
            assert key in expected
            assert getattr(association, key) == expected[key]


def test_species_parse_gene(panther_species_genes):
    for gene_info, expected in panther_species_genes:
        species, gene_id = parse_gene_info(gene_info, 
                                           panther_taxon_map, 
                                           db_to_curie_map)
                                           ### NEEDS FIXING / MANUAL CURATION (TODO STILL)
                                           ###make_ncbi_taxon_gene_map("./data/gene_info.gz", relevant_ncbi_cols, relevant_ncbi_taxons))
        
        # Assert that the parsed species and gene ID match the expected values
        assert species in panther_taxon_map
        assert gene_id in expected # Allows for muliple mapping values to be present
        

def test_exclude_species_parse_gene(exclude_species_genes):
    for gene_info in exclude_species_genes:
        species, gene_id = parse_gene_info(gene_info, 
                                           panther_taxon_map, 
                                           db_to_curie_map)

        # Assert that the species and gene ID are empty for excluded species
        assert not species
        assert not gene_id
