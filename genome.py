"""
Methods for retrieving and interacting with whole genome data.

All methods use NCBI Entrez and GenBank file formats.
"""

from typing import List
from xml.etree import ElementTree

from Bio import SeqIO, Entrez


"""
Find a given gene in genome annotation tree

:param genome: Genomic `SeqIO` object

:param term: name or id of gene to search for in annotation tree

:returns: `SeqIO.features` data for matching value
"""
def find_in_genome(genome, term):
    for feature in genome.features:
        if 'gene' in feature.qualifiers.keys() and term in feature.qualifiers['gene']:
            return feature


"""
Get RefSeq genomes from NCBI using a given query.

:param organism: search query to find organism on GenBank

:returns: List of `SeqRecord` objects containing reference genomes
"""
def get_genome_list(organism: str) -> List['SeqIO']:
    # TODO: have more expansive and comprehensive query
    query = '"%s"[Organism] AND (Refseq[Filter] AND "bioproject nuccore"[Filter] AND "scope monoisolate"[Filter])' \
            'NOT plasmid[filter]' % organism
    with Entrez.esearch(db='bioproject', term=query, retmode='text') as handle:
        records = Entrez.read(handle)

    accessions = []
    for uid in records['IdList']:
        with Entrez.efetch(db='bioproject', id=uid) as handle:
            root = ElementTree.parse(handle).getroot()
            accession = root[0][0][0][0].attrib['accession']
            accessions.append(accession)

    ids = []
    for accession in accessions:
        with Entrez.esearch(db='nuccore', term=accession) as handle:
            ids.extend(Entrez.read(handle)['IdList'])

    genomes = []
    for i in ids:
        genomes.append(fetch_genome(i))

    # Exclude any sequence that has "plasmid" in title
    return [i for i in genomes if 'plasmid' not in i.description]


"""
Lookup and download genome by `id`.

This returns a list of `SeqIO` in the event that multiple genomes are returned, but they *should* not

:param uid: NCBI id of genome to retrieve

:returns: `SeqIO` corresponding to matching genome
"""
def fetch_genome(uid: str) -> SeqIO:
    with Entrez.efetch(db='nuccore', id=uid, rettype='gbwithparts', retmode='text') as handle:
        return SeqIO.read(handle, 'gb')


"""
Search for genomes using accession number, then retrieve using id.

This returns a list of `SeqIO` in the event that multiple genomes are returned, but they *should* not.

:param accession: Genome search term to use

:returns list of `SeqIO` objects that *should* correspond to given accession number
"""
def fetch_genome_by_accession(accession: str) -> List['SeqIO']:

    with Entrez.esearch(db='nuccore', term=accession) as handle:
        ids = Entrez.read(handle)['IdList']

    genomes = []
    for i in ids:
        genomes.append(fetch_genome(i))
    return genomes