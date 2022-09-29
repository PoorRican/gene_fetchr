# fetcher
Library to search and retrieve genomes, genes, vectors and other genomic data from public data sources.
Currently, Genbank is the only supported data source. Scrapers for Addgene are planned for future releases.

A list of all genes can be extracted from a genome sequence annotation feature.
Using the ASN.1 file format, more detailed information (i.e: description, summary) than what the `gb` file format
offers can be extracted. In future releases, data such as gene ontology can be extracted.

Generated data is flat and unprocessed and can be used to convert to other data formats, has applications for
computational biology, or be fed to ML algorithms.