"""Compare the species order baked into generated kernels with a transliteration
of the C library's Markowitz ordering (sol-max src/reorder.c, astra-xhigh
src/linear.c), for configurations run through copies of the generators.

Usage: kernel_order_check.py ENTRY GENERATED_C CONFIG_JSON
Reads the configuration with a minimal CAMP reader (JSON only, as in this run).
"""
import json, re, sys
from pathlib import Path

MASK = (1 << 64) - 1

def camp(config):
    config = Path(config)
    root = json.loads(config.read_text())
    data = list(root.get('camp-data', []))
    for name in root.get('camp-files', []):
        data += json.loads((config.parent / name).read_text()).get('camp-data', [])
    species = [d['name'] for d in data if d['type'] == 'CHEM_SPEC']
    reactions = []
    for m in (d for d in data if d['type'] == 'MECHANISM'):
        for r in m['reactions']:
            t = r['type']
            if t == 'SURFACE':
                rs, ps = [r['gas-phase reactant']], list(r.get('gas-phase products', {}))
            elif t == 'EMISSION':
                rs, ps = [], [r['species']]
            elif t == 'FIRST_ORDER_LOSS':
                rs, ps = [r['species']], []
            else:
                rs = []
                for n, o in r.get('reactants', {}).items():
                    rs += [n] * int(o.get('qty', 1))
                ps = list(r.get('products', {}))
            if t == 'BRANCHED':
                reactions.append((rs, list(r['alkoxy products'])))
                reactions.append((rs, list(r['nitrate products'])))
            else:
                reactions.append((rs, ps))
    return species, reactions

def c_sol_max(species, reactions):
    n = len(species); v = {s: i for i, s in enumerate(species)}
    p = [[0] * n for _ in range(n)]
    for i in range(n): p[i][i] = 1
    for rs, ps in reactions:
        for a in rs:
            for b in rs + ps: p[v[b]][v[a]] = 1
    alive = [1] * n; perm = []
    for _ in range(n):
        best, cost = n, None
        for x in range(n):
            if not alive[x]: continue
            rows = sum(1 for j in range(n) if alive[j] and p[x][j])
            cols = sum(1 for j in range(n) if alive[j] and p[j][x])
            cand = ((rows - 1) * (cols - 1)) & MASK
            if best == n or cand < cost: best, cost = x, cand
        perm.append(best); alive[best] = 0
        for row in range(n):
            if alive[row] and p[row][best]:
                for col in range(n):
                    if alive[col] and p[best][col]: p[row][col] = 1
    return [species[i] for i in perm]

def c_astra_xhigh(species, reactions):
    n = len(species); v = {s: i for i, s in enumerate(species)}
    g = [[0] * n for _ in range(n)]
    for rs, ps in reactions:
        for a in rs:
            for b in rs + ps: g[v[b]][v[a]] = 1
    row = [sum(g[i]) for i in range(n)]; col = [sum(g[i][j] for i in range(n)) for j in range(n)]
    alive = [1] * n; perm = []
    for _ in range(n):
        pivot, cost = None, None
        for i in range(n):
            if alive[i]:
                cand = ((row[i] - 1) * (col[i] - 1)) & MASK
                if pivot is None or cand < cost: pivot, cost = i, cand
        perm.append(pivot); alive[pivot] = 0
        for i in range(n):
            if alive[i] and g[i][pivot]:
                for j in range(n):
                    if alive[j] and g[pivot][j] and not g[i][j]:
                        g[i][j] = 1; row[i] += 1; col[j] += 1
        for i in range(n):
            if alive[i]:
                if g[pivot][i]: col[i] -= 1
                if g[i][pivot]: row[i] -= 1
    return [species[i] for i in perm]

def kernel_order(entry, text, species, reactions):
    """Recover species -> variable index from the generated forcing kernel."""
    if entry == 'sol-max':
        body = text.split('_force_0(')[1]
        blocks = re.findall(r'double rate = k\[(\d+)\*s\+c\];(.*?)\n  \}', body, re.S)
        pat = r'rate \*= y\[(\d+)\*s\+c\]'
    else:
        body = text.split('static void forcing(')[1].split('static void jacobian(')[0]
        blocks = re.findall(r'double r=rates\[(\d+)\*cells\+c\];(.*?)\n    \}', body, re.S)
        pat = r'r\*=y\[(\d+)\*cells\+c\]'
    index = {}
    for k, block in blocks:
        rs = reactions[int(k)][0] if entry == 'astra-xhigh' else None
        ys = [int(x) for x in re.findall(pat, block)]
        if rs is not None:
            for name, y in zip(rs, ys):
                if index.setdefault(name, y) != y: raise SystemExit(f'inconsistent index for {name}')
    return index, blocks, pat

entry, generated, config = sys.argv[1:4]
species, reactions = camp(config)
text = Path(generated).read_text()
c_order = (c_sol_max if entry == 'sol-max' else c_astra_xhigh)(species, reactions)
if entry == 'sol-max':
    # sol-max sorts reactions by rate type before emitting kernels; recover the
    # emitted order from the loader's chem_model_reaction calls.
    calls = re.findall(r'const size_t reactants\[\] = \{([^}]*)\}', text)
    emitted = [[species[int(x)] for x in c.split(',') if x.strip()] for c in calls]
    nreact = [int(x) for x in re.findall(r'chem_model_reaction\(m, reactants, (\d+),', text)]
    emitted = [e[:n] for e, n in zip(emitted, nreact)]
    body = text.split('_force_0(')[1]
    blocks = re.findall(r'double rate = k\[(\d+)\*s\+c\];(.*?)\n  \}', body, re.S)
    index = {}
    for k, block in blocks:
        ys = [int(x) for x in re.findall(r'rate \*= y\[(\d+)\*s\+c\]', block)]
        for name, y in zip(emitted[int(k)], ys):
            if index.setdefault(name, y) != y: raise SystemExit(f'inconsistent index for {name}')
else:
    index, _, _ = kernel_order(entry, text, species, reactions)
kernel = [None] * len(species)
for name, y in index.items(): kernel[y] = name
print(f'{entry} {Path(config).parent.name}: kernel order {kernel}')
print(f'{entry} {Path(config).parent.name}: C order      {c_order}')
known = [k for k in kernel if k is not None]
print('AGREE' if all(kernel[i] in (None, c_order[i]) for i in range(len(c_order))) else 'DIFFER',
      f'({len(known)} of {len(species)} species appear as reactants)')
