import collections
import csv
import json
from pathlib import Path
import sys
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT/'tools'))
from import_ts1 import translate
from config import read_mechanism
from generate import generate

class Ts1Translation(unittest.TestCase):
    def test_translation_and_reference_defaults(self):
        source = (ROOT/'../micm/test/benchmark/ts1_mechanism.hpp').read_text()
        species, reactions, defaults = translate(source)
        configured_species, configured_reactions = read_mechanism(ROOT/'configs/ts1/config.json')
        config = json.loads((ROOT/'configs/ts1/config.json').read_text())
        self.assertEqual((species, reactions, defaults),
                         (configured_species, configured_reactions, config['custom parameter defaults']))
        self.assertEqual(collections.Counter(r['type'] for r in reactions),
                         {'ARRHENIUS': 361, 'USER_DEFINED': 142, 'TROE': 31, 'SURFACE': 13})
        with (ROOT/'../reference/ts1_defaults.csv').open(newline='') as stream:
            reference = next(csv.DictReader(stream))
        actual = {s['name']: s.get('initial concentration [mol m-3]', 0)
                  for s in species if not s.get('third body')}
        actual.update(defaults)
        actual['temperature'] = config['temperature [K]']
        actual['pressure'] = config['pressure [Pa]']
        self.assertEqual(set(actual), set(reference))
        for name, value in reference.items():
            with self.subTest(name=name):
                self.assertEqual(actual[name], float(value))
        output = generate(ROOT/'configs/ts1/config.json')
        self.assertIn('const ChemMechanism generated_mechanism = {209, 547,', output)

if __name__ == '__main__':
    unittest.main()
