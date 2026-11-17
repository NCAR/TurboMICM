import pathlib
import sys
import unittest

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parents[1] / 'tools'))
from config import read_mechanism
from generate import generate

ROOT = pathlib.Path(__file__).resolve().parents[2] / 'micm/test/tutorial/configs'

class ConfigurationTests(unittest.TestCase):
    def test_all_examples_read(self):
        for name, extension, n_species, n_reactions in [
            ('robertson', 'json', 3, 3),
            ('rate_constants_no_user_defined', 'yaml', 7, 6),
            ('rate_constants_user_defined', 'json', 7, 9),
        ]:
            with self.subTest(name=name):
                species, reactions = read_mechanism(ROOT / name / ('config.' + extension))
                self.assertEqual((len(species), len(reactions)), (n_species, n_reactions))

    def test_generator_rejects_unsupported_rate(self):
        with self.assertRaisesRegex(ValueError, 'ARRHENIUS'):
            generate(ROOT / 'rate_constants_no_user_defined/config.yaml')

if __name__ == '__main__':
    unittest.main()
