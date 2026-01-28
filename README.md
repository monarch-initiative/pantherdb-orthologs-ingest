# pantherdb-orthologs-ingest

Koza ingest for PantherDB ortholog data, transforming gene orthology relationships into Biolink model format.

## Data Source

[PantherDB](http://www.pantherdb.org/) provides gene orthology data derived from phylogenetic trees. The AllOrthologs file contains pairwise ortholog relationships between genes across multiple species.

Data is downloaded from: `http://data.pantherdb.org/ftp/ortholog/current_release/AllOrthologs.tar.gz`

Additionally, NCBI gene_info is used for identifier mapping: `https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_info.gz`

## Output

This ingest produces:
- **GeneToGeneHomologyAssociation** - Orthology relationships between genes from different species

## Supported Species

The ingest processes orthology data for these species:
- Human, Mouse, Rat
- Zebrafish, Chicken, Xenopus
- Dog, Cow, Pig
- Fly, Worm, Yeast
- Dictyostelium, Schizosaccharomyces pombe, Aspergillus

## Usage

```bash
# Install dependencies
just install

# Run full pipeline
just run

# Or run steps individually
just download      # Download PantherDB data
just transform-all # Run Koza transform
just test          # Run tests
```

## Requirements

- Python 3.10+
- [uv](https://github.com/astral-sh/uv) package manager
- [just](https://github.com/casey/just) command runner

## License

MIT
