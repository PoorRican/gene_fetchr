"""
Methods to retrieve files from NCBI Genbank using Entrez
"""

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
def get_genome_list(organism: str):
    """
    Get RefSeq genomes from NCBI using a given query.
    :param organism: search query

    :returns List of `SeqRecord` objects containing reference genome
    """
    # TODO: have more expansive query
    query = '"%s"[Organism] AND (Refseq[Filter] AND "bioproject nuccore"[Filter] AND "scope monoisolate"[Filter])' \
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
    for id in ids:
        print("Trying %s" % id)
        with Entrez.efetch(db='nuccore', id=id, rettype='gb', retmode='text') as handle:
            genomes.append(SeqIO.read(handle, 'gb'))

    # Exclude any sequence that has "plasmid" in title
    return [i for i in genomes if 'plasmid' not in i.description]


# Get genome data by accession number
def fetch_genome(id):
    return NotImplemented

# Fetch all features found within a genome
def get_features(genome):
    return NotImplemented

