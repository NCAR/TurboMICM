"""Generate a test registry, rejecting duplicate or misspelled reference IDs."""
from pathlib import Path
import re
import sys

known={line.split()[0] for line in Path('tests/cases.txt').read_text().splitlines() if line and not line.startswith('#')}
cases={}
for path in sorted(Path('tests').glob('test_*.c')):
    for prefix,suite,name in re.findall(r'\bTEST\((\w+),\s*(\w+),\s*(\w+)\)',path.read_text()):
        identifier=f'{prefix}:{suite}.{name}'
        if identifier not in known or identifier in cases:
            raise ValueError('Unknown or duplicate test: '+identifier)
        cases[identifier]=f'test_{prefix}_{suite}_{name}'
    for identifier in re.findall(r'\bPROGRAM\((\w+)\)',path.read_text()):
        if identifier not in known or identifier in cases:
            raise ValueError('Unknown or duplicate program: '+identifier)
        cases[identifier]=f'test_program_{identifier}'
lines=['/* Generated registry of implemented cases. */']
for function in cases.values(): lines.append(f'extern void {function}(void);')
lines.append('static const struct test_case cases[]={')
for identifier,function in cases.items(): lines.append(f'    {{"{identifier}", {function}, NULL}},')
lines.extend(['};',''])
Path(sys.argv[1]).parent.mkdir(parents=True,exist_ok=True)
Path(sys.argv[1]).write_text('\n'.join(lines))
