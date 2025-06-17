"""
Utility functions for Panther Orthology data processing
"""
import os
import gzip
from collections import Counter



# These names should match pantherdb shorthand names for each species
# Example... https://www.pantherdb.org/genomes/genome.jsp?taxonId=9606 --> Short Name: HUMAN
panther_taxon_map = {"HUMAN": "9606",
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
                     "YEAST": "4932"}
                     # Additional species for future here...
                     # "FELCA": "9685",  # Felis catus - domestic cat

relevant_ncbi_cols = ['#tax_id', 
                      'GeneID', 
                      'Symbol', 
                      'LocusTag', 
                      'Synonyms', 
                      'dbXrefs', 
                      'Symbol_from_nomenclature_authority',
                      'Full_name_from_nomenclature_authority',
                      'Other_designations']

# Entries with Gene/Orthology identifier namespaces that need modifying to match our CURIEs
# Keys are the pantherdb namespace, and values are the CURIE namespace 
# (Note, many key/value pairs are the same for simplicty in downstream processing
db_to_curie_map = {"HGNC":"HGNC",
                   "MGI":"MGI",
                   "RGD":"RGD",
                   "SGD":"SGD",
                   "ZFIN":"ZFIN",
                   "dictyBase":"dictyBase",
                   "PomBase":"PomBase",
                   "Xenbase":"Xenbase",
                   
                   "FlyBase":"FB",      
                   "WormBase":"WB",     
                   "Ensembl": "ENSEMBL"}
                   #"EnsemblGenome": "ENSEMBL",    # TO DO: Needs reviewing still if we need to map to ncbi or ensembl

# Used in make_ncbi_taxon_gene_map function to filter for only species we are interested in
relevant_ncbi_taxons = {v:'' for v in panther_taxon_map.values()} 



def make_ncbi_taxon_gene_map(gene_info_file: str, relevant_columns: list, taxon_catalog: dict):
    
    # Ensure relevant columns has #tx_id as the first entry
    if relevant_columns[0] != "#tax_id":
        raise RuntimeError("- '#tax_id' must be first element present in relevant_columns arg... Exiting")
    
    # We don't want entries equivalant to this from this file
    exclude_terms = {"-":''}
    
    # Don't want to use pandas here (for memory and other reasons relating to speed)
    taxa_gene_map = {tx_id:{} for tx_id in taxon_catalog}   # Many-->1 mapping dictionary
    taxa_remove_map = {tx_id:Counter() for tx_id in taxon_catalog} # Removes unreliable mapping keys from taxa_gene_map
    
    
    with gzip.open(gene_info_file, 'rt') as infile:
        
        # Read header line into memory to index relevant column fields
        hinfo = {hfield:i for i,hfield in enumerate(infile.readline().strip('\r').strip('\n').split('\t'))}
        
        ccc = 0
        # Now loop through each line, and create a map back to taxon / NCBIGene:xyz ...
        for line in infile:
            cols = line.strip('\r').strip('\n').split('\t')
            rel_data = [str(cols[hinfo[r]]) for r in relevant_columns]
            tx_id = rel_data[0]
            ncbi_gene_id = cols[hinfo["GeneID"]]
            
            # Only consume species we are interested in
            if tx_id not in taxon_catalog:
                continue
            
            # Find reliable mapping keys to this NCBI gene id
            # We take the set() of relevant mapping keys here... that if the same_id is reported on the same line
            # This removes the possibility of removing an id that is reported twice on the same line
            removed = []
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
                        taxa_gene_map[tx_id].update({key_to_ncbi:ncbi_gene_id})
                    else:
                        taxa_remove_map[tx_id][key_to_ncbi] += 1
            
    
    # Remove unreliable mapping keys to this NCBI gene id
    for tx_id in taxa_remove_map:
        for remove_key, rcount in taxa_remove_map[tx_id].items():
            ##print("- Taxon {} | Removing {} | Count {}".format(tx_id, remove_key, rcount))
            del taxa_gene_map[tx_id][remove_key]
    
    # Return cleaned map back to a ncbi gene id that can be normalized later down road
    return taxa_gene_map


def parse_gene_info(gene_info, taxon_map, curie_map, fallback_map):
    """
    This function takes a panther gene information string and returns the species name and gene identifier
    in a standardized format by converting to CURIEs based on predefined mappings, and using the 
    uniprotkb id as a fallback. We also remove ensemble version/transcript ids from the tail end of ensembl ids,
    and we also filter out species that are not in our taxon map. Below are examples of the transformation process

    HUMAN|HGNC=16435|UniProtKB=Q8TCT9 --> HUMAN, HGNC:16435
    SCHPO|PomBase=SPBP23A10.09|UniProtKB=Q9P7X6 --> SCHPO, PomBase:SPBP23A10.09
    CANLF|Ensembl=ENSCAFG00845009646.1|UniProtKB=A0A8I3N1X7 --> CANLF, Ensembl:ENSCAFG00845009646
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
    matched = 0
    fback = 0
    unikb = 0
    
    # Check if gene id can be mapped directly to kg build preffered gene ids
    if gene_split[0] in curie_map:
        gene = "{}:{}".format(curie_map[gene_split[0]], gene_split[-1]) # We use -1 here to avoid things like MGI=MGI=95886
        matched = 1
        
    # Otherwise, fall back onto ncbi gene map if possible, 
    elif gene_split[1] in fallback_map[taxon_id]:
        g_id = fallback_map[taxon_id][gene_split[1]]
        gene = "NCBIGene:{}".format(g_id)
        fback = 1
    
    # Use uniprotkb id as last resort and format e.g. UniProtKB=Q8TCT9 => UniProtKB:Q8TCT9
    else:
        gene = "{}".format(cols[-1].replace("=", ":"))
        unikb += 1

        
    # Lastly we need to strip version numbers off from ENSEMBL IDs,
    # (e.g. ENSG00000123456.1 => ENSG00000123456)
    if gene.startswith("ENSEMBL:") and (":ENS" in gene):
        gene = gene.split(".")[0]
    
    return species, gene,