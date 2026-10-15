"""Reader validation independent of the numerical tutorial equivalence cases."""
import json
import pathlib
import tempfile
import unittest

import configure


class ReaderTests(unittest.TestCase):
    def read(self, data):
        with tempfile.TemporaryDirectory(dir='build') as directory:
            path = pathlib.Path(directory) / 'config.json'
            path.write_text(json.dumps(data))
            return configure.read(path)

    def minimal(self):
        return {'concentration units': 'mol m-3', 'camp-data': [
            {'type': 'CHEM_SPEC', 'name': 'A'}, {'type': 'CHEM_SPEC', 'name': 'B'},
            {'type': 'MECHANISM', 'name': 'test', 'reactions': [
                {'type': 'ARRHENIUS', 'reactants': {'A': {}}, 'products': {'B': {}}}]}]}

    def test_tutorials(self):
        for directory, extension in [('robertson','json'), ('rate_constants_user_defined','json'),
                                      ('rate_constants_no_user_defined','yaml')]:
            mechanism = configure.read(pathlib.Path('mechanisms/tutorial') / directory / ('config.'+extension))
            self.assertTrue(mechanism.reactions)
            self.assertIn('const ck_mechanism', configure.emit(mechanism, 'fixture'))

    def test_duplicate_species(self):
        data = self.minimal(); data['camp-data'][1]['name']='A'
        with self.assertRaises(configure.ConfigurationError): self.read(data)

    def test_unknown_species(self):
        data = self.minimal(); data['camp-data'][2]['reactions'][0]['reactants']={'missing':{}}
        with self.assertRaises(configure.ConfigurationError): self.read(data)

    def test_nonfinite(self):
        data = self.minimal(); data['camp-data'][2]['reactions'][0]['A']=float('inf')
        with self.assertRaises(configure.ConfigurationError): self.read(data)

    def test_qty(self):
        for value in [True, 0, -1, 1.5]:
            data = self.minimal(); data['camp-data'][2]['reactions'][0]['reactants']={'A':{'qty':value}}
            with self.assertRaises(configure.ConfigurationError): self.read(data)

    def test_units_and_duplicates(self):
        data = self.minimal(); data['concentration units']='molecule cm-3'
        reaction=data['camp-data'][2]['reactions'][0]
        reaction.update(reactants=['A','A'], products=[{'species':'B','yield':2}, {'species':'B','yield':.5}], A=2)
        mechanism=self.read(data); parsed=mechanism.reactions[0]
        self.assertEqual(parsed['reactants'],[0,0]); self.assertEqual(parsed['products'],[(1,2),(1,.5)])
        self.assertEqual(parsed['rate']['A'],2*configure.AVOGADRO_CM)

    def test_extensions(self):
        for kind, constant in [('TERNARY_CHEMICAL_ACTIVATION_JPL19','CK_TERNARY_JPL19'),
                               ('TAYLOR_SERIES','CK_TAYLOR'),('REVERSIBLE','CK_REVERSIBLE')]:
            data=self.minimal(); data['camp-data'][2]['reactions'][0]['type']=kind
            mechanism=self.read(data); self.assertIn(constant,configure.emit(mechanism,'fixture'))
            if kind.endswith('JPL19'): self.assertEqual(mechanism.reactions[0]['rate']['k0_D'],298)
            if kind=='TAYLOR_SERIES': self.assertEqual(mechanism.reactions[0]['rate']['coefficients'],[1])

    def test_cycle(self):
        with tempfile.TemporaryDirectory(dir='build') as directory:
            path=pathlib.Path(directory)/'config.json';path.write_text(json.dumps({'camp-files':['config.json']}))
            with self.assertRaises(configure.ConfigurationError): configure.read(path)

    def test_bad_symbol(self):
        mechanism=self.read(self.minimal())
        with self.assertRaises(configure.ConfigurationError): configure.emit(mechanism,'invalid-name')

    def test_missing_surface_properties(self):
        data=self.minimal(); data['camp-data'][2]['reactions'][0].update(type='SURFACE',label='surface',**{'gas-phase reactant':'A'})
        with self.assertRaises(configure.ConfigurationError): self.read(data)


if __name__ == '__main__':
    unittest.main()
