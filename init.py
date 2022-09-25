from Bio import Entrez

FILENAME = 'email'


def set_email() -> None:
    """
    Sets email value for `Entrez` module functionality.

    This should be called near the top of any module that uses `Entrez`
    """
    with open(FILENAME, 'r') as handle:
        Entrez.email = handle.read()
