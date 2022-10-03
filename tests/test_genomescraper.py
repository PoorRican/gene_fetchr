import unittest
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from genbank.GenomeScraper import GenomeScraper

FILE = 'sequences/U00096.gb'


class TestGenomeScraper(unittest.TestCase):
    def test_functor(self):
        """
        Test that class correctly acts as a functor
        """
        with open(FILE, 'r') as f:
            gb = SeqIO.read(f, 'gb')
            genome = GenomeScraper(gb)

        for attr in ('keys', 'values', 'features', 'scrape'):
            self.assertTrue(hasattr(genome, attr), '`GenomeScraper` doesn\'t have \'%s\' attribute' % attr)

    def test_invalid_features(self):
        """
        Test that error is raised when less than 2 features occur
        """
        seq = 'ACTCTCTCTCTCTCT'
        genome = SeqRecord(Seq(seq))
        with self.assertRaises(ValueError):
            GenomeScraper(genome)

    def test_properties(self):
        """
        Test that all property functions return correct data
        """
        with open(FILE, 'r') as f:
            gb = SeqIO.read(f, 'gb')
            genome = GenomeScraper(gb)

        # test string values
        for attr in ('topology', 'molecule_type', 'organism', 'comment', 'id', 'name'):
            self.assertTrue(hasattr(genome, attr))
            self.assertIsInstance(getattr(genome, attr), str)

    def test_update_gene(self):
        """
        Test that `_update_gene` works correctly by passing unique data,
        then checking to see that data was not overwritten.
        """
        return NotImplemented
