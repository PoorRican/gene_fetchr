"""
Methods to retrieve files from NCBI Genbank using Entrez
"""

from typing import List
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
def get_genome_list(organism: str) -> List[SeqIO]:
    """
    Get RefSeq genomes from NCBI using a given query.
    :param organism: search query

    :returns List of `SeqRecord` objects containing reference genome
    """
    # TODO: have more expansive query
    query = '"%s"[Organism] AND (Refseq[Filter] AND "bioproject nuccore"[Filter] AND "scope monoisolate"[Filter]) NOT plasmid[filter]' \
            % organism
    with Entrez.esearch(db='bioproject', term=query, retmode='text') as handle:
        records = Entrez.read(handle)

    accessions = []
    for id in records['IdList']:
        with Entrez.efetch(db='bioproject', id=id) as handle:
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
def fetch_genome(id: str) -> SeqIO:
    with Entrez.efetch(db='nuccore', id=id, rettype='gb', retmode='text') as handle:
        return SeqIO.read(handle, 'gb')

def fetch_gene(id: str) -> SeqIO:
    with Entrez.efetch(db='gene', id=id, rettype='gbwithparts', retmode='text') as handle:
        return SeqIO.read(handle, 'gb')


# Get genome data by accession number
def fetch_genome_by_accession(accession: str) -> List[SeqIO]:
    """
    Search for genomes using accession number, then retrieve using id
    :param accession: Genome search term to use

    :returns list of `SeqIO` objects
    """
    with Entrez.esearch(db='nuccore', term=accession) as handle:
        ids = Entrez.read(handle)['IdList']

    genomes = []
    for i in ids:
        genomes.append(fetch_genome(i))
    return genomes


# Fetch all features found within a genome

def get_all_genes(genome: SeqIO) -> List[SeqIO]:
    # raise error if there are no annotations
    if len(genome.features) <= 1:
        raise ValueError('Genome does not have any gene annotations')

    # check if first feature is `source`
    start = 0
    if (genome.features[0].type == 'source'): start = 1

    genes = []
    # fetch all genes
    for gene in genome.features:
        gid = [i[i.index(':')+1:] for i in genomes[0].features[1].qualifiers['db_xref'] if 'GeneID' in i][0]
        gene = fetch_gene(gid)
        genes.append(gene)
    return genes

def get_gene_metadata(id: str) -> dict['title': str, 'summary': str]:
    """
    Get non-essential metadata not found in GenBank files by using defailt ASN.1 specification.

    :param id: GeneID
    :return: dict containing `summary` and `title`
    """
    with Entrez.efetch(db='gene', id=id, retmode='xml') as handle:
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