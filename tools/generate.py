"""Generate a C mass-action mechanism from a CAMP configuration.

Currently only photolysis reactions are supported. Rate constants are supplied
by the caller in configuration order. Unsupported types fail explicitly.
"""
import argparse
import json
from config import read_mechanism


def generate(path):
    species, reactions = read_mechanism(path)
    indices = {s["name"]: i for i, s in enumerate(species)}
    reactants, products, rows = [], [], []
    for reaction in reactions:
        if reaction["type"] != "PHOTOLYSIS":
            raise ValueError(f"unsupported reaction type: {reaction['type']}")
        rb, pb = len(reactants), len(products)
        for name, properties in reaction.get("reactants", {}).items():
            qty = properties.get("qty", 1)
            if not isinstance(qty, int) or qty < 1:
                raise ValueError("reactant qty must be a positive integer")
            reactants.extend([indices[name]] * qty)
        for name, properties in reaction.get("products", {}).items():
            products.append((indices[name], float(properties.get("yield", 1))))
        rows.append((rb, len(reactants)-rb, pb, len(products)-pb))
    lines = ['#include "chem.h"', '']
    lines.append('const char *const generated_species[] = {' +
                 ', '.join(json.dumps(s["name"], ensure_ascii=True) for s in species) + '};')
    lines.append('static const size_t reactants[] = {' +
                 ', '.join(map(str, reactants or [0])) + '};')
    lines.append('static const ChemProduct products[] = {' +
                 ', '.join('{%d, %.17g}' % p for p in (products or [(0, 0)])) + '};')
    lines.append('static const ChemReaction reactions[] = {' +
                 ', '.join('{%d, %d, %d, %d}' % row for row in (rows or [(0, 0, 0, 0)])) + '};')
    lines.append('const ChemMechanism generated_mechanism = {%d, %d, reactions, reactants, products};' %
                 (len(species), len(reactions)))
    return '\n'.join(lines) + '\n'


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('config')
    parser.add_argument('output')
    args = parser.parse_args()
    with open(args.output, 'w', encoding='utf-8') as out:
        out.write(generate(args.config))
