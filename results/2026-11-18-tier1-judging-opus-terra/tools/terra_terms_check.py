"""Map the species indices in a Terra entry's generated C back to names and
compare the reactions, as multisets, with the configuration's.

    terra_terms_check.py ENTRY GENERATED_C CONFIG_JSON
Third bodies and parameterized species are dropped from the configuration side,
since the generated tables carry them as rate multipliers.
"""
import json, re, sys
from collections import Counter
from pathlib import Path

def num(s):
    s = s.strip()
    return float.fromhex(s) if 'x' in s.lower() else float(s)

def parse(entry, t):
    if entry == 'terra-max':
        names = re.findall(r'^\s*\{ "([^"]*)", ', re.search(r'ChemSpecies species\[\] = \{(.*?)\n\};', t, re.S).group(1), re.M)
        R = [int(x) for x in re.findall(r'-?\d+', re.search(r'reactants\[\] = \{([^}]*)\}', t).group(1))]
        P = [(int(i), num(y)) for i, y in re.findall(r'\{ (\d+), ([^ }]+) \}', re.search(r'ChemProduct products\[\] = \{(.*?)\};', t).group(1))]
        rows = [tuple(map(int, m)) for m in re.findall(r'^\s*\{ (\d+), (\d+), (\d+), (\d+), \d+, \{', t, re.M)]
    elif entry == 'terra-xhigh':
        names = re.findall(r'"([^"]*)"', re.search(r'species_names\[\] = \{(.*?)\};', t, re.S).group(1))
        R = [int(i) for i, _ in re.findall(r'\{(-?\d+), ([^}]+)\}', re.search(r'ChemTerm reactants\[\] = \{(.*?)\n\};', t, re.S).group(1))]
        P = [(int(i), num(y)) for i, y in re.findall(r'\{(-?\d+), ([^}]+)\}', re.search(r'ChemTerm products\[\] = \{(.*?)\n\};', t, re.S).group(1))]
        rows = [tuple(map(int, m)) for m in re.findall(r'\{CHEM_\w+, (\d+), (\d+), (\d+), (\d+), -?\d+, \{', t)]
    else:
        names = re.findall(r'"([^"]*)"', re.search(r'generated_species\[\] = \{([^}]*)\}', t).group(1))
        R = [int(x) for x in re.findall(r'\d+', re.search(r'size_t reactants\[\] = \{([^}]*)\}', t).group(1))]
        P = [(int(i), num(y)) for i, y in re.findall(r'\{(\d+), ([^}]+)\}', re.search(r'ChemProduct products\[\] = \{(.*?)\};', t).group(1))]
        rows = [tuple(map(int, m)) for m in re.findall(r'\{(\d+), (\d+), (\d+), (\d+)\}', re.search(r'ChemReaction reactions\[\] = \{(.*?)\};', t).group(1))]
    out = []
    for rb, rn, pb, pn in rows:
        rs = [names[i] for i in R[rb:rb+rn] if i >= 0]
        ps = [(names[i], round(y, 12)) for i, y in P[pb:pb+pn]]
        out.append((tuple(sorted(rs)), tuple(sorted(ps))))
    return names, out

def expected(config):
    root = Path(config)
    cfg = json.loads(root.read_text())
    species, reactions = [], []
    for f in cfg['camp-files']:
        for e in json.loads((root.parent / f).read_text())['camp-data']:
            if e['type'] == 'CHEM_SPEC': species.append(e)
            elif e['type'] == 'MECHANISM': reactions += e['reactions']
    drop = {s['name'] for s in species if s.get('third body') or 'parameterization' in s}
    out = []
    for r in reactions:
        typ = r['type']
        if typ == 'BRANCHED':
            forms = [(r['reactants'], r['alkoxy products']), (r['reactants'], r['nitrate products'])]
        elif typ == 'SURFACE':
            forms = [({r['gas-phase reactant']: {}}, r.get('gas-phase products', {}))]
        elif typ == 'FIRST_ORDER_LOSS':
            forms = [({r['species']: {}}, {})]
        elif typ == 'EMISSION':
            forms = [({}, {r['species']: {}})]
        else:
            forms = [(r.get('reactants', {}), r.get('products', {}))]
        for rs, ps in forms:
            a = sorted(n for n, p in rs.items() if n not in drop for _ in range(p.get('qty', 1)))
            b = sorted((n, round(float(p.get('yield', 1)), 12)) for n, p in ps.items() if n not in drop)
            out.append((tuple(a), tuple(b)))
    return out

entry, path, config = sys.argv[1:4]
names, got = parse(entry, Path(path).read_text())
want = expected(config)
ok = Counter(got) == Counter(want)
print(f"{entry} {Path(path).name}: {len(got)} reactions, species {names}: {'match' if ok else 'MISMATCH'}")
if not ok:
    print('  generated only:', list((Counter(got) - Counter(want)).elements()))
    print('  configuration only:', list((Counter(want) - Counter(got)).elements()))
