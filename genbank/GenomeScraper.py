"""
Methods to retrieve and interact with gene data.

`get_gene_metadata` uses NCBI ASN.1 BioSeq specification since GenBank file format
returns no metadata in their "gene" db.
"""

from collections.abc import Generator
from Bio import SeqRecord, SeqFeature
from collections import UserDict

from GeneScraper import GeneScraper


class GenomeScraper(UserDict):
    """
    Functor that accepts a `SeqRecord` object and aggregates and fetches metadata of all `SeqFeature` objects.

    Todo:
        - Separate queries to enable multithreading
    """

    def __init__(self, genome: SeqRecord):
        """

        Args:
            genome: `SeqRecord` of interest
                `SeqRecord` must have annotations (other than `source`), otherwise
        """

        super().__init__()

        if len(genome.features) <= 1:
            raise ValueError('Genome does not have any gene annotations')

        self.genome = genome
        self.scrape()

    @staticmethod
    def _get_gid(gene: SeqFeature) -> str:
        """
        Get `SeqRecord` using `db_xref`.

        If `GeneID` is not present, use all values in `db_xref` as search terms to get `GeneID`

        Args:
            gene: feature to fetch from NCBI

        Returns:
            valid GeneID to pass to `Entrez.efetch`

        Todo:
            This should be moved to `GeneScraper`
        """

        gid = [i[i.index(':') + 1:] for i in gene.qualifiers['db_xref'] if 'GeneID' in i][0]
        return gid

    @property
    def features(self) -> Generator['SeqRecord']:
        """
        Iterator for `self.genome`

        This function skips the `source` feature, since it is unnecessary.

        Yields:
            Members of `self.genome.features`
        """
        # check if first feature is `source`
        if self.genome.features[0].type == 'source':
            start = 1
        else:
            start = 0

        for gene in self.genome.features[start:]:
            yield gene

    def _update_gene(self, gene: SeqFeature) -> None:
        """
        Update `self.data` with `gene`.

        This function is necessary since genomes have duplicate, explicit `SeqFeatures` for
        gene and CDS data.

        Args:
            gene:

        Returns:

        """
        return NotImplemented

    def scrape(self) -> None:
        """
        Populate `self.data` with `GeneScraper` objects for each feature in genome
        """
        for gene in self.features:
            data = GeneScraper(gene)
            # TODO: update `self.data` instead of overwriting
            self.data[data.id] = data
