"""
Methods to retrieve and interact with gene data.

Both `fetch_gene` and `get_all_genes` use Entrez to actively download GenBank data onto memory,
and use the GenBank file format.

Conversely, `get_gene_metadata` uses NCBI's ASN.1 BioSeq specification since more information is stored in this type.
"""

from typing import List
from Bio import Entrez, SeqIO
from xml.etree import ElementTree


def fetch_gene(uid: str) -> SeqIO:
    with Entrez.efetch(db='gene', id=uid, rettype='gbwithparts', retmode='text') as handle:
        return SeqIO.read(handle, 'gb')


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

    :paramater uid: GeneID

    :returns: dict containing `summary` and `title`
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
