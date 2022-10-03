from collections.abc import Generator
from Bio import SeqRecord, SeqFeature, Seq
from collections import UserDict

from GeneScraper import GeneScraper


class GenomeScraper(UserDict):
    """
    Functor that scrapes all useful data on all genes within an annotated genome.
    Gene data is presented via a dict-like interface.

    A passthrough for underlying `SeqRecord` properties are also provided.

    Todo:
        - Separate queries to enable multithreading
    """

    def __init__(self, genome: SeqRecord):
        """
        Populate instance `data` with features from `genome`.

        Args:
            genome: `SeqRecord` of interest
                `SeqRecord` must have annotations (other than `source`), otherwise `ValueError` is raised.
        """

        super().__init__()

        if len(genome.features) <= 1:
            raise ValueError('Genome does not have any gene annotations')

        self.genome = genome
        self.scrape()

    @property
    def topology(self) -> str:
        return self.genome.annotations['topology']

    @property
    def molecule_type(self) -> str:
        return self.genome.annotations['molecule_type']

    @property
    def organism(self) -> str:
        return self.genome.annotations['organism']

    @property
    def comment(self) -> str:
        return self.genome.annotations['comment']

    @property
    def id(self) -> str:
        return self.genome.id

    @property
    def name(self):
        return self.genome.name

    @property
    def sequence(self) -> Seq:
        return self.genome.seq

    @property
    def description(self) -> str:
        return self.genome.description

    def __repr__(self):
        if self.description:
            return self.description
        elif self.id:
            return self.id
        elif self.name:
            return self.name
        elif self.topology and self.molecule_type:
            return '%s %s with %d sub-features' % (self.topology, self.molecule_type, len(self.data))

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
            data = GeneScraper(gene, self.genome)
            # TODO: update `self.data` instead of overwriting
            self.data[data.id] = data
