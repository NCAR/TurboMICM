"""Map the species indices in each generated file back to names.

Usage: terms_check.py ENTRY GENERATED_C
Prints the reactions as the generated C names them, with third bodies dropped
where the entry drops them, so they can be compared with the configuration.
"""
import re, sys
from pathlib import Path

def ints(s):
    return [int(x) for x in re.findall(r'-?\d+', s)]

def pairs(s):
    return [int(a) for a, _ in re.findall(r'\{\s*(\d+)u?\s*,\s*([^}]*)\}', s)]

def parse(entry, t):
    out = []
    if entry == 'sol-high':
        names = re.findall(r'"([^"]*)"', re.search(r'species_names\[\] = \{([^}]*)\}', t).group(1))
        for i in range(len(re.findall(r'reactants_\d+\[\]', t))):
            r = ints(re.search(rf'reactants_{i}\[\] = \{{([^}}]*)\}}', t).group(1))
            p = pairs(re.search(rf'products_{i}\[\] = \{{(.*?)\}};', t).group(1))
            n = re.search(r'\{(\d+),(\d+),reactants_%d,' % i, t)
            out.append((r[:int(n.group(1))], p[:int(n.group(2))]))
    elif entry == 'sol-xhigh':
        names = re.findall(r'\.name="([^"]*)", \.molecular', t)
        for line in re.findall(r'\{\.nreactants=.*', t):
            nr, rv, np_, pv = re.search(r'nreactants=(\d+), \.reactants=(\S+), \.nproducts=(\d+), \.products=(\S+),', line).groups()
            r = ints(re.search(re.escape(rv) + r'\[\] = \{([^}]*)\}', t).group(1)) if rv != 'NULL' else []
            p = pairs(re.search(re.escape(pv) + r'\[\] = \{(.*?)\};', t).group(1)) if pv != 'NULL' else []
            out.append((r, p))
    elif entry == 'sol-max':
        names = re.findall(r'chem_model_species\(m, "([^"]*)"', t)
        for rs, ps, nr, np_ in re.findall(r'reactants\[\] = \{([^}]*)\}, products\[\] = \{([^}]*)\};.*?chem_model_reaction\(m, reactants, (\d+), products, yields, (\d+),', t, re.S):
            out.append((ints(rs)[:int(nr)], ints(ps)[:int(np_)]))
    elif entry == 'astra-high':
        names = re.findall(r'\{\.name="([^"]*)"', t)
        n = len(re.findall(r'\.rate = chem_rate_default', t))
        for i in range(n):
            m = re.search(rf'reactants_{i}\[\] = \{{([^}}]*)\}}', t)
            q = re.search(rf'products_{i}\[\] = \{{(.*?)\}};', t)
            out.append((ints(m.group(1)) if m else [], pairs(q.group(1)) if q else []))
    elif entry == 'astra-xhigh':
        names = re.findall(r'\{\.name="([^"]*)", \.type', t)
        for nr, np_, i in re.findall(r'\{(\d+),(\d+),r(\d+),p\d+,', t):
            r = ints(re.search(rf'r{i}\[\]=\{{([^}}]*)\}}', t).group(1))[:int(nr)]
            p = pairs(re.search(rf'p{i}\[\]=\{{(.*?)\}};', t).group(1))[:int(np_)]
            out.append((r, p))
    elif entry == 'astra-max':
        names = re.findall(r'ch_model_add_species\(m,"([^"]*)"', t)
        for nr, rs, np_, ps in re.findall(r'ch_model_add_reaction\(m,(\d+),(\(int\[\]\)\{[^}]*\}|NULL),(\d+),(\(int\[\]\)\{[^}]*\}|NULL),', t):
            out.append((ints(rs.replace('(int[])', '')) if rs != 'NULL' else [], ints(ps.replace('(int[])', '')) if ps != 'NULL' else []))
    elif entry == 'luna-max':
        names = re.findall(r'\{\.name="([^"]*)", \.molecular', t)
        R = ints(re.search(r'chem_reactant_data\[\] = \{([^}]*)\}', t).group(1))
        P = ints(re.search(r'chem_product_data\[\] = \{([^}]*)\}', t).group(1))
        for ro, rc, po, pc in re.findall(r'\.reactant_offset=(\d+), \.reactant_count=(\d+), \.product_offset=(\d+), \.product_count=(\d+)', t):
            out.append((R[int(ro):int(ro)+int(rc)], P[int(po):int(po)+int(pc)]))
    elif entry == 'luna-xhigh':
        names = re.findall(r'\{ "([^"]*)", [^,]*, [^,]*, [^,]*, [^,]*, 0u \}', t)
        R = pairs(re.search(r'generated_reactants\[\d+\] = \{(.*?)\n\};', t, re.S).group(1))
        P = pairs(re.search(r'generated_products\[\d+\] = \{(.*?)\n\};', t, re.S).group(1))
        for ro, rc, po, pc in re.findall(r'\{ (\d+)u, (\d+)u, (\d+)u, (\d+)u, \{ \.kind', t):
            out.append((R[int(ro):int(ro)+int(rc)], P[int(po):int(po)+int(pc)]))
    elif entry == 'luna-high':
        names = re.findall(r'\{"([^"]*)", \d+, ', t)
        for kind, rp, nr, pp, np_ in re.findall(r'\{(CHEM_\w+), "[^"]*", (\S+), (\d+), (\S+), (\d+),', t):
            r = pairs(re.search(re.escape(rp) + r'\[\] = \{(.*?)\};', t, re.S).group(1)) if rp != 'NULL' else []
            p = pairs(re.search(re.escape(pp) + r'\[\] = \{(.*?)\};', t, re.S).group(1)) if pp != 'NULL' else []
            out.append((r, p))
    return names, out

entry, path = sys.argv[1:3]
names, reactions = parse(entry, Path(path).read_text())
print(f'{entry}: species table {names}')
for k, (r, p) in enumerate(reactions):
    try:
        print(f'  R{k}: ' + ' + '.join(names[i] for i in r) + ' -> ' + ' + '.join(names[i] for i in p))
    except IndexError as e:
        print(f'  R{k}: index out of range: reactants {r} products {p}')
