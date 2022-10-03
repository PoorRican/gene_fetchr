import unittest
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
from genbank.GeneScraper import GeneScraper


class TestGeneScraper(unittest.TestCase):
    def test_sequence(self):
        seq = 'ACTCTCTCTCTCT'
        genome = SeqRecord(Seq(seq))

        # test top strand
        feature = SeqFeature(FeatureLocation(0, 5, 1), strand=1, type='gene')
        scraped = GeneScraper(feature, genome)
        self.assertEqual(scraped.sequence, genome[0:5].seq)

        # test complimentary strand
        feature.strand = -1
        scraped = GeneScraper(feature, genome)
        self.assertEqual(scraped.sequence, genome[0:5].reverse_complement().seq)

    def test_minimal_data(self):
        """Test that there are minimal requirements for `GeneScraper`"""
        empty = SeqFeature()
        no_strand = SeqFeature(FeatureLocation(0, 5), type='test')
        no_type = SeqFeature(FeatureLocation(0, 5), strand=1)

        with self.assertRaises(ValueError):
            GeneScraper(empty)
            GeneScraper(no_strand)
            GeneScraper(no_type)
