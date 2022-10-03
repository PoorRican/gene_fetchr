import unittest
from ontologies.OntologyScraper import OntologyScraper


mak = "P23917"


class TestOntologyScraper(unittest.TestCase):
    def test_class_attributes(self):
        """Test that `u` and `go` are class attributes"""
        self.assertTrue(hasattr(OntologyScraper, 'u'))
        self.assertTrue(hasattr(OntologyScraper, 'go'))

    def test_functor_functionality(self):
        """Test that class can be called"""
        o = OntologyScraper()
        self.assertIsNotNone(o(mak), 'OntologyScraper.__call__ did not work as expected')

    def test_get_ontologies(self):
        """Verify values return from `get_ontologies`"""
        o = OntologyScraper()
        ontologies = o(mak)
        self.assertIsInstance(ontologies, list, 'get_ontologies did not return list')

        keys = ('name', 'id', 'definition', 'aspect')
        for ont in ontologies:
            for key in keys:
                self.assertTrue(key in ont.keys())
                self.assertIsNotNone(ont[key])


if __name__ == '__main__':
    unittest.main()
