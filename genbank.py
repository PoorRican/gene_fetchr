"""
Methods to retrieve files from NCBI Genbank using Entrez
"""

from typing import List, TypeVar
from Bio import Entrez, SeqIO
from xml.etree import ElementTree


# Iterate through `seq.features`
def find_in_genome(genome, term):
    """
    Find a given gene in genome
    :param genome: `SeqIO` object
    :param term: gene to find
    :return:
    """
    for feature in genome.features:
        if 'gene' in feature.qualifiers.keys() and term in feature.qualifiers['gene']:
            return feature


# Get a list of reference genomes matching organism
def get_genome_list(organism: str) -> List['SeqIO']:
    """
    Get RefSeq genomes from NCBI using a given query.

    :param organism: search query to find organism on GenBank

    :returns: List of `SeqRecord` objects containing reference genome
    """
    # TODO: have more expansive query
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


# Get genome data by id
def fetch_genome(uid: str) -> SeqIO:
    with Entrez.efetch(db='nuccore', id=uid, rettype='gb', retmode='text') as handle:
        return SeqIO.read(handle, 'gb')


def fetch_gene(uid: str) -> SeqIO:
    with Entrez.efetch(db='gene', id=uid, rettype='gbwithparts', retmode='text') as handle:
        return SeqIO.read(handle, 'gb')


# Get genome data by accession number
def fetch_genome_by_accession(accession: str) -> List['SeqIO']:
    """
    Search for genomes using accession number, then retrieve using id.

    This returns a list of `SeqIO` in the event that multiple genomes are returned, but they *should* not

    :param accession: Genome search term to use

    :returns list of `SeqIO` objects
    """
    with Entrez.esearch(db='nuccore', term=accession) as handle:
        ids = Entrez.read(handle)['IdList']

    genomes = []
    for i in ids:
        genomes.append(fetch_genome(i))
    return genomes


def get_all_genes(genome: SeqIO) -> List['SeqIO']:
    """
    Fetch all annotated genes in a genome from GenBank

    :param genome: Annotated genomic `SeqIO` of interest

    :return: List of annotated sequences
    """
    # raise error if there are no annotations
    if len(genome.features) <= 1:
        raise ValueError('Genome does not have any gene annotations')

    # check if first feature is `source`
    if genome.features[0].type == 'source':
        start = 1
    else:
        start = 0

    genes = []
    # fetch all genes
    for gene in genome.features[start:]:
        gid = [i[i.index(':') + 1:] for i in gene.qualifiers['db_xref'] if 'GeneID' in i][0]
        data = fetch_gene(gid)
        genes.append(data)
    return genes


def get_gene_metadata(uid: str) -> dict[str, str]:
    """
    Get non-essential metadata not found in GenBank files by using defailt ASN.1 specification.

    :param uid: GeneID
    :return: dict containing `summary` and `title`
    """
    with Entrez.efetch(db='gene', id=uid, retmode='xml') as handle:
        parsed = ElementTree.parse(handle)
    root = parsed.getroot()

    try:
        prot_desc = root[0].find('Entrezgene_prot').find('Prot-ref').find('Prot-ref_desc').text
    except ValueError:
        prot_desc = None

    try:
        summary = root[0].find('Entrezgene_summary').text
    except ValueError:
        summary = None

    return {
        'summary': summary,
        'prot_desc': prot_desc,
    }
