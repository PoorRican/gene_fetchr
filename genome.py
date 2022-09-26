"""
Methods for retrieving and interacting with whole genome data.

All methods use NCBI Entrez and GenBank file formats.
For more details on genes or products, refer to `gene.py` or `protein.py`
"""

from typing import List
from xml.etree import ElementTree
from warnings import warn

from Bio import SeqIO, Entrez, SeqRecord, SeqFeature


def find_in_genome(genome, term) -> SeqFeature:
    """
    Find a given gene in genome annotation tree

    :param genome: Genomic `SeqRecord`

    :param term: name or id of gene to search for in annotation tree

    :returns: `SeqFeature` data. Contains GeneID and location
    """
    for feature in genome.features:
        if 'gene' in feature.qualifiers.keys() and term in feature.qualifiers['gene']:
            return feature


def get_genome_list(organism: str) -> List['SeqRecord']:
    """
    Get RefSeq genomes from NCBI using a given query.

    :param organism: search query to find organism on GenBank

    :returns: List of reference genome `SeqRecord` objects
    """

    # TODO: have more expansive and comprehensive query
    query = '"%s"[Organism] AND (Refseq[Filter] AND "bioproject nuccore"[Filter] AND "scope monoisolate"[Filter])' \
            'NOT plasmid[filter]' % organism
    with Entrez.esearch(db='bioproject', term=query, rettype='gb', retmode='text') as handle:
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

    return genomes


def fetch_genome(gid: str) -> SeqRecord:
    """
    Lookup and download genome by `gi`.

    This returns a list of `SeqIO` in the event that multiple genomes are returned, but they *should* not

    :param gid: NCBI id (`gi`) of genome to retrieve. From `Entrez.esearch` or ASN.1 data...

    :returns: `SeqIO` corresponding to matching genome
    """
    with Entrez.efetch(db='nuccore', id=gid, rettype='gbwithparts', retmode='text') as handle:
        return SeqIO.read(handle, 'gb')


def fetch_genome_by_accession(accession: str) -> SeqRecord:
    """
    Search for genomes using accession number, then retrieve using `gi`

    :raises A warning is raised if accession number returned 0 or more than 1 `SeqFeatures`.
    The first returned `GeneID` is used in the event that more than one record is returned by `esearch`.

    :param accession: Genome search term to use

    :returns `SeqRecord` matching accession number
    """

    with Entrez.esearch(db='nuccore', term=accession) as handle:
        uid = Entrez.read(handle)['IdList']

    if len(uid) == 0:
        raise ValueError("Invalid genome accession number (%s)")
    if len(uid) != 1:
        warn("Genome accession number (%s) returned more than one GeneID value" % accession)
        return fetch_genome(uid[0])
    else:
        return fetch_genome(uid)
