# fetcher
Library to search and retrieve genomes, genes, vectors and other genomic data.

Has knowledge of data-type and can translate accordingly (ie: gene from genome, protein from gene, etc.).
Feeds data to `constructor` and `builder`.

Uses NCBI Entrez and Addgene as data sources, but in reality this module is simply hardcoded retrieval methods to
pivot Genbank.