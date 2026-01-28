"""
Utility functions for Panther Orthology data processing
"""

# These names should match pantherdb shorthand names for each species
# Example... https://www.pantherdb.org/genomes/genome.jsp?taxonId=9606 --> Short Name: HUMAN
panther_taxon_map = {
    "HUMAN": "9606",
    "MOUSE": "10090",
    "CANLF": "9615",   # Canis lupus familiaris - domestic dog
    "BOVIN": "9913",   # Bos taurus - cow
    "PIG": "9823",     # Sus scrofa - pig
    "RAT": "10116",
    "CHICK": "9031",
    "XENTR": "8364",   # Xenopus tropicalis - tropical clawed frog
    "DANRE": "7955",
    "DROME": "7227",
    "CAEEL": "6239",
    "DICDI": "44689",
    "EMENI": "227321",  # Emericella nidulans (strain FGSC A4 etc.) (Aspergillus nidulans)
    "SCHPO": "4896",
    "YEAST": "4932",
    # Additional species for future here...
    # "FELCA": "9685",  # Felis catus - domestic cat
}

# Entries with Gene/Orthology identifier namespaces that need modifying to match our CURIEs
# Keys are the pantherdb namespace, and values are the CURIE namespace
# (Note, many key/value pairs are the same for simplicity in downstream processing)
db_to_curie_map = {
    "HGNC": "HGNC",
    "MGI": "MGI",
    "RGD": "RGD",
    "SGD": "SGD",
    "ZFIN": "ZFIN",
    "dictyBase": "dictyBase",
    "PomBase": "PomBase",
    "Xenbase": "Xenbase",
    "FlyBase": "FB",
    "WormBase": "WB",
    "Ensembl": "ENSEMBL",
    # Note: EnsemblGenome identifiers are handled by falling back to NCBI gene map
    # or UniProtKB, since they're in the symbol namespace not ENSEMBL gene namespace
}


def parse_gene_info(gene_info, taxon_map, curie_map, fallback_map=None, koza_transform=None):
    """
    This function takes a panther gene information string and returns the species name and gene identifier
    in a standardized format by converting to CURIEs based on predefined mappings, and using the
    uniprotkb id as a fallback. We also remove ensemble version/transcript ids from the tail end of ensembl ids,
    and we also filter out species that are not in our taxon map. Below are examples of the transformation process

    HUMAN|HGNC=16435|UniProtKB=Q8TCT9 --> HUMAN, HGNC:16435
    SCHPO|PomBase=SPBP23A10.09|UniProtKB=Q9P7X6 --> SCHPO, PomBase:SPBP23A10.09
    CANLF|Ensembl=ENSCAFG00845009646.1|UniProtKB=A0A8I3N1X7 --> CANLF, Ensembl:ENSCAFG00845009646

    Args:
        gene_info: Panther gene information string (species|gene|uniprotkb_id)
        taxon_map: Dictionary mapping panther species names to NCBI taxon IDs
        curie_map: Dictionary mapping panther database prefixes to CURIE prefixes
        fallback_map: Nested dict {taxon_id: {identifier: ncbi_gene_id}} for fallback lookups (legacy)
        koza_transform: KozaTransform object for lookups via koza mappings (preferred)
    """

    cols = gene_info.split("|") # species|gene|uniprotkb_id
    species = cols[0]

    # Exit condition (saves compute when there are many rows to process..)
    if species not in taxon_map:
        return None, None

    # Now assign our gene to its "rightful" prefix... If no reasonable prefix exists (HGNC, MGI, etc.),
    # then we use the UniprotKB ID prefix as a fallback. Connections can be rescued through
    # normalization process via uniprotkb protein ids

    # Our prefered order, is Organism specific (HGNC, PomBase, ZFIN)
    taxon_id = taxon_map[species]
    gene_split = cols[1].split("=")

    # Check if gene id can be mapped directly to kg build preffered gene ids
    if gene_split[0] in curie_map:
        gene = "{}:{}".format(curie_map[gene_split[0]], gene_split[-1]) # We use -1 here to avoid things like MGI=MGI=95886

    # Otherwise, fall back onto ncbi gene map if possible
    else:
        # Try koza lookup first (production mode)
        if koza_transform is not None:
            composite_key = f"{taxon_id}|{gene_split[1]}"
            ncbi_gene_id = koza_transform.lookup(composite_key, "ncbi_gene_id", map_name="ncbi_gene")
            if ncbi_gene_id != composite_key:  # lookup returns key if not found
                gene = ncbi_gene_id
            else:
                # Use uniprotkb id as last resort
                gene = "{}".format(cols[-1].replace("=", ":"))
        # Legacy fallback_map mode (for tests)
        elif fallback_map is not None and gene_split[1] in fallback_map.get(taxon_id, {}):
            g_id = fallback_map[taxon_id][gene_split[1]]
            gene = "NCBIGene:{}".format(g_id)
        # Use uniprotkb id as last resort
        else:
            gene = "{}".format(cols[-1].replace("=", ":"))


    # Lastly we need to strip version numbers off from ENSEMBL IDs,
    # (e.g. ENSG00000123456.1 => ENSG00000123456)
    if gene.startswith("ENSEMBL:") and (":ENS" in gene):
        gene = gene.split(".")[0]

    return species, gene