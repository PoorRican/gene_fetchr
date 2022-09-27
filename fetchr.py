from Bio import SeqIO
from GenomeScraper import GenomeScraper


FILE = 'sequences/U00096.gb'

if __name__ == '__main__':
    with open(FILE, 'r') as file:
        gb = SeqIO.read(file, 'gb')
        genome = GenomeScraper(gb)