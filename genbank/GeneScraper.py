from Bio import SeqFeature, Entrez, SeqRecord, Seq
from collections import UserDict
from xml.etree import ElementTree
from typing import Union, Sized


def _deconstruct(item: Union[list, set, tuple]) -> Union[object, list, set, tuple]:
    """
    Deconstruct iterable containers into single value if there is only one member of list,
    otherwise, preserve iterable container.

    Meant to break out values from containers when there is only one object. It cannot be assumed
    that parsed files will always have one object in container for any 1 datatype. In `SeqFeature`,
    containerization of a single object seems to only occur in the `qualifiers` container.

    Args:
        item: Any iterable container (i.e: tuple, list, set)

    Returns:
        sole member of container, or the container is untouched
    """

    if len(item) == 1:
        return item[0]
    else:
        return item


class GeneScraper(UserDict):
    def __init__(self, gene: SeqFeature, genome: SeqRecord = None):
        super().__init__()

        # store for debugging
        self.gene = gene
        self.genome = genome

        if gene.type == 'gene':
            self.data = self._from_gene(gene)
        elif gene.type == 'CDS':
            self.data = self._from_cds(gene)
        elif gene.type == 'ncRNA':
            self.data = self._from_nc_rna(gene)
        elif gene.type == 'mobile_element':
            self.data = self._from_mobile_element(gene)
        else:
            self.data = self._from_gene(gene)

    def __repr__(self):
        return self.id

    @property
    def id(self) -> Union[str, int]:
        if 'gene' in self.data.keys():
            return self.data['gene']
        elif 'name' in self.data.keys():
            return self.data['name']
        elif self.data['type'] == 'mobile_element':
            return self.data['mobile_element_type']
        else:
            return "Unknown feature (%d, %d)" % (self.data['start'], self.data['end'])

    @property
    def sequence(self) -> Seq:
        """Get nucleotide sequence from parent genome.

        This functionality is deferred to `biopython` for simplicity.

        Returns
            `Seq` string defined at this gene's location.
        """
        return self.gene.extract(self.genome.seq)

    @staticmethod
    def _qualifiers(gene: SeqFeature) -> dict:
        """
        Get generic `qualifiers` data

        Keys include 'locus_tag', 'db_xref', 'name', 'product', and 'note'.
        If any key is not found, it is skipped.
        """
        keys = ('gene', 'gene_synonym', 'locus_tag', 'db_xref', 'name', 'product', 'note', 'mol_type')
        data = {}

        for key in keys:
            if key in gene.qualifiers:
                data[key] = _deconstruct(gene.qualifiers[key])

        return data

    @staticmethod
    def _minimal_data(gene: SeqFeature) -> dict:
        """
        Ensure that minimal data requirements are met.
        """
        if gene.location is None or gene.type is None or gene.location.strand is None or \
           gene.location.start is None or gene.location.end is None:
            raise ValueError('Minimal feature attributes not met')
        data = {
            'strand': gene.location.strand,
            'type': gene.type,
            'start': int(gene.location.start),
            'end': int(gene.location.end),
        }
        return data

    @staticmethod
    def _from_gene(gene: SeqFeature) -> dict:
        """
        Extract data from a feature.

        This is called by default for any `SeqFeature` not caught by `__init__`
        (i.e: 'misc_feature', 'tRNA', 'rRNA').

        Args:
            gene: feature to extract data from

        Returns:
            dict containing keys pulled from `SeqFeature`.
            Serves as transient datatype before being added to `self.data`
        """

        if gene:            # get rid of static check warning
            data = {}
            data.update(GeneScraper._minimal_data(gene))
            data.update(GeneScraper._qualifiers(gene))
            return data

    @staticmethod
    def _from_cds(gene: SeqFeature) -> dict:
        """
        Extract data from `cds` feature type

        Args:
            gene: feature to extract data from

        Returns:
            `CdsData` pulled from `SeqFeature`.
            Serves as transient datatype before being added to `self.data`
        """

        # TODO: in what case would there be more than one codon_start/codon_end/product?
        keys = ('codon_start', 'transl_table', 'protein_id', 'translation')
        data = {}

        for key in keys:
            if key in gene.qualifiers:
                data[key] = _deconstruct(gene.qualifiers[key])

        data.update(GeneScraper._from_gene(gene))
        return data

    @staticmethod
    def _from_nc_rna(gene: SeqFeature) -> dict:
        data = {
            'ncRNA_class': _deconstruct(gene.qualifiers['ncRNA_class']),
        }

        data.update(GeneScraper._from_gene(gene))
        return data

    @staticmethod
    def _from_mobile_element(gene: SeqFeature) -> dict:
        data = {'mobile_element_type': _deconstruct(gene.qualifiers['mobile_element_type'])}

        data.update(GeneScraper._from_gene(gene))
        return data

    @staticmethod
    def _gene_metadata(gid: str) -> dict[str, str]:
        """
        Get non-essential metadata not found in GenBank files by using ASN.1 spec.

        Since querying NCBI is time-consuming, this function is not called by `__init__`

        Args:
            gid: GeneID

        Returns:
            dict containing `summary` and `title`

        Todo:
            Return a namedtuple
        """
        with Entrez.efetch(db='gene', id=gid, retmode='xml') as handle:
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
