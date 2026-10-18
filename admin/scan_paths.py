"""List commands in a Codex session record that touch paths outside the run's own inputs.

    scan_paths.py ROLLOUT RUN_DIR
"""
import json, re, sys
path, run_dir = sys.argv[1], sys.argv[2].rstrip('/')
root = run_dir.rsplit('/', 1)[0]
allowed = (run_dir, root + '/micm', root + '/reference', root + '/scenarios')
hits = []
for line in open(path):
    d = json.loads(line)
    p = d.get('payload') if isinstance(d.get('payload'), dict) else {}
    if p.get('type') != 'item_completed':
        continue
    item = p.get('item', {})
    if item.get('type') != 'CommandExecution':
        continue
    cmd = ' '.join(item.get('command', []))
    for m in re.finditer(r'(/Users/fillmore[^\s\'";|&)]*|\.\./[^\s\'";|&)]*)', cmd):
        target = m.group(1)
        if target.startswith('../'):
            name = target[3:].split('/')[0]
            if name in ('micm', 'reference', 'scenarios', ''):
                continue
        elif target.startswith(allowed) or target.startswith(('/Users/fillmore/miniforge3', '/Users/fillmore/.codex')):
            continue
        hits.append((d['timestamp'][11:19], target, cmd[-160:].replace('\n', ' ')))
for h in hits:
    print(*h, sep='  ')
print(f'{len(hits)} references outside the allowed inputs')
