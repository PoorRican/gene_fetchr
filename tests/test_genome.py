import unittest

from gene_fetchr import set_email
from genbank.genome import *


set_email()


class TestGenomeModule(unittest.TestCase):
    def test_find_in_feature(self):
        term = "mak"
        with open('sequences/U00096.gb', 'r') as file:
            genome = SeqIO.read(file, 'gb')
        feature = find_in_genome(genome, term)
        self.assertIn('mak', feature.qualifiers['gene'], "Incorrect gene found or improper gene format")

    def test_get_genome_list(self):
        """ Test that `get_genome_list` gets appropriate reference genome
        """
        accession = "NC_000913"     # reference genome of E. coli str. K-12 substr. MG1655 as of 9/24/2022

        genomes = get_genome_list('e coli')
        accessions = [i.name for i in genomes]
        self.assertIn(accession, accessions, "Query does not get correct assertion")

    def test_fetch_genome(self):
        gid = "556503834"           # found via `Entrez.esearch` or from `seq-set > seq > id > gi` in ASN.1 data

        accession = "NC_000913"
        genome = fetch_genome(gid)
        self.assertEqual(accession, genome.name, "Incorrect record when fetching by id")

    def test_fetch_genome_by_accession(self):
        accession = "NC_000913.3"
        name = "NC_000913"

        genome = fetch_genome_by_accession(accession)
        self.assertEqual(genome.name, name, "Incorrect record when fetching by accession number")

        # test that error is thrown with incorrect accession number
        with self.assertRaises(ValueError):
            incorrect = 'aoeuaoeu_incorrect accession number. no way esearch will return data...'
            fetch_genome_by_accession(incorrect)

        # TODO: test that warning is raised when more than one `SeqRecord` is returned by `esearch`
        # https://docs.python.org/3/library/warnings.html#testing-warnings


if __name__ == '__main__':
    unittest.main()
