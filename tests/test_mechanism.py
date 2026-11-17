import csv
import json
import tempfile
import unittest
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools'))
from mechanism import MechanismError, read_mechanism
from generate import compile_mechanism

ROOT = Path(__file__).resolve().parents[1]
TUTORIAL = ROOT / 'tests/tutorial_configs'


class ReaderTests(unittest.TestCase):
    def test_tutorials(self):
        yaml_mech = read_mechanism(TUTORIAL / 'rate_constants_no_user_defined/config.yaml')
        json_mech = read_mechanism(TUTORIAL / 'rate_constants_user_defined/config.json')
        robertson = read_mechanism(TUTORIAL / 'robertson/config.json')
        self.assertEqual((len(yaml_mech['species']), len(yaml_mech['reactions'])), (7, 6))
        self.assertEqual((len(json_mech['species']), len(json_mech['reactions'])), (7, 9))
        self.assertEqual((len(robertson['species']), len(robertson['reactions'])), (3, 3))
        self.assertEqual(yaml_mech['reactions'], json_mech['reactions'][:6])

    def test_tutorial_rate_unit_conversion(self):
        mech = read_mechanism(TUTORIAL / 'rate_constants_no_user_defined/config.yaml')
        reactions = compile_mechanism(mech)[-1]
        factor = 1e-6 * 6.02214076e23
        ternary = next(r for r in reactions if r[0] == 'TERNARY')[-1]
        troe = next(r for r in reactions if r[0] == 'TROE')[-1]
        self.assertEqual(ternary[0], 1.2)
        self.assertEqual(ternary[3], 2.6 / factor)
        self.assertEqual(troe[0], 1.8e-30 * factor**2)
        self.assertEqual(troe[3], 2.8e-11 * factor)

    def test_ts1_counts_and_defaults(self):
        mech = read_mechanism(ROOT / 'configs/ts1/config.json')
        self.assertEqual((len(mech['species']), len(mech['reactions'])), (210, 547))
        self.assertEqual(sum(bool(s.get('third body')) for s in mech['species']), 1)
        self.assertEqual(mech['initial conditions']['jo2_a'], 1)
        self.assertEqual(mech['initial conditions']['O2'], 8.90261411)

    def test_ts1_defaults_match_reference(self):
        mech = read_mechanism(ROOT / 'configs/ts1/config.json')
        with (ROOT / 'tests/ts1_defaults.csv').open(newline='') as stream:
            reference = next(csv.DictReader(stream))
        defaults = mech['initial conditions']
        for name, value in reference.items():
            if name in ('temperature', 'pressure'):
                self.assertEqual(mech['initial environment'][name], float(value), name)
                continue
            self.assertEqual(defaults.get(name, 0.0), float(value), name)

    def test_parameterized_species_are_not_state_variables(self):
        mech = {'species': [{'type': 'CHEM_SPEC', 'name': 'A'},
                            {'type': 'CHEM_SPEC', 'name': 'M', 'parameterization': {'air density': 2}}],
                'reactions': [{'type': 'ARRHENIUS', 'reactants': {'M': {}}, 'products': {'A': {}}}],
                'initial conditions': {}, 'rate units': 'mol m-3'}
        species, _, reactants, _, _ = compile_mechanism(mech)
        self.assertEqual([s['name'] for s in species], ['A'])
        self.assertEqual(reactants, [(-2, 1.0)])

    def test_bad_species_reference(self):
        with tempfile.TemporaryDirectory(dir=ROOT / 'build') as directory:
            path = Path(directory)
            (path / 'config.json').write_text(json.dumps({'camp-files': ['data.json']}))
            (path / 'data.json').write_text(json.dumps({'camp-data': [
                {'type': 'CHEM_SPEC', 'name': 'A'},
                {'type': 'MECHANISM', 'reactions': [{'type': 'ARRHENIUS', 'reactants': {'B': {}}}]}]}))
            with self.assertRaises(MechanismError):
                read_mechanism(path / 'config.json')

    def test_bad_initial_environment(self):
        with tempfile.TemporaryDirectory(dir=ROOT / 'build') as directory:
            path = Path(directory)
            (path / 'data.json').write_text(json.dumps({'camp-data': [{'type': 'CHEM_SPEC', 'name': 'A'}]}))
            for environment in ({'temperature': -1}, {'pressure': True}, {'unknown': 1}, []):
                (path / 'config.json').write_text(json.dumps({'camp-files': ['data.json'], 'initial environment': environment}))
                with self.assertRaises(MechanismError):
                    read_mechanism(path / 'config.json')


if __name__ == '__main__':
    unittest.main()
