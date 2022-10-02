from bioservices import UniProt, QuickGO
from typing import List


class OntologyScraper(object):
    """
    Singleton that retrieves ontologies for a given protein product.
    Once this class is instantiated, the class method `get_ontologies` is used to fetch ontology data.

    Data is retrieved by cross-referencing `UniProtKB` data with `QuickGO`,
    therefore ontologies can only be found for coding sequences. These queries also mean that
    fetching is very slow, so any functionality that incorporates returned data should
    be asynchronous.

    Notes
        This class uses `bioservices`, which is licensed under GPL which requires that source code
        be publicly available.
    """

    # store these on a class level because they are computationally expensive to instantiate
    u = UniProt(verbose=False)
    go = QuickGO()

    @classmethod
    def __call__(cls, uniprot_id: str) -> List[dict]:
        return cls.get_ontologies(uniprot_id)

    @classmethod
    def get_ontologies(cls, uniprot_id: str) -> List[dict]:
        """
        Return GO terms for a given protein.

        First, UniProtKB is queried, getting GO URI's. Then details of those URIs
        are cross-referenced using QuickGO. This process is very slow, and must be
        implemented asynchronously.

        Args:
            uniprot_id: UniProt id for given product

        Returns:
            List containing retrieved ontologies for given gene

        """
        results = cls.u.retrieve(uniprot_id)
        go_ids = []
        for i in results['uniProtKBCrossReferences']:
            if i['database'] == 'GO':
                go_ids.append(i['id'])

        ontologies = cls.go.get_go_terms(','.join(go_ids))
        results = []
        for i in ontologies:
            results.append({
                'name': i['name'],
                'id': i['id'],
                'definition': i['definition'],
                'aspect': i['aspect'],
            })
        return results
