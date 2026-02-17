# PantherDB Orthologs

[PantherDB](http://www.pantherdb.org/) (Protein Analysis Through Evolutionary Relationships) provides gene orthology data derived from phylogenetic trees. The AllOrthologs file contains pairwise ortholog relationships between genes across multiple species.

Data is downloaded from: `http://data.pantherdb.org/ftp/ortholog/current_release/AllOrthologs.tar.gz`

Additionally, [NCBI gene_info](https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz) is used for identifier mapping when direct database mappings are not available.

### Species Included

- Human (NCBITaxon:9606)
- Mouse (NCBITaxon:10090)
- Rat (NCBITaxon:10116)
- Dog (NCBITaxon:9615)
- Cow (NCBITaxon:9913)
- Pig (NCBITaxon:9823)
- Chicken (NCBITaxon:9031)
- Frog (NCBITaxon:8364)
- Zebrafish (NCBITaxon:7955)
- Fruit fly (NCBITaxon:7227)
- C. elegans (NCBITaxon:6239)
- Dictyostelium (NCBITaxon:44689)
- Aspergillus nidulans (NCBITaxon:227321)
- Fission yeast (NCBITaxon:4896)
- Baker's yeast (NCBITaxon:4932)

## Ortholog Associations

Gene identifiers from the PantherDB file are mapped to standard CURIE namespaces (HGNC, MGI, RGD, SGD, ZFIN, FB, WB, PomBase, dictyBase, Xenbase, ENSEMBL). When a direct mapping is not available, NCBI gene IDs are used as a fallback via a lookup table, with UniProtKB as a last resort. ENSEMBL version numbers are stripped (e.g. `ENSG00000123456.1` becomes `ENSG00000123456`).

Rows where either species is not in the supported taxon set are skipped.

**Biolink Captured:**

- `biolink:GeneToGeneHomologyAssociation`
    - id (UUID)
    - subject (gene CURIE)
    - predicate (`biolink:orthologous_to`)
    - object (gene CURIE)
    - has_evidence (PANTHER.FAMILY ID)
    - aggregator_knowledge_source (`["infores:monarchinitiative"]`)
    - primary_knowledge_source (`infores:panther`)
    - knowledge_level (`knowledge_assertion`)
    - agent_type (`not_provided`)

## Citation

Thomas PD, Ebert D, Muruganujan A, Mushayahama T, Albou LP, Mi H. PANTHER: Making genome-scale phylogenetics accessible to all. Protein Science. 2022;31(1):8-22. doi: 10.1002/pro.4218. PMID: 34717010

## License

MIT
