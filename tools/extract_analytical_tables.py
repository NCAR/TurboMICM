"""Retain the numerical reference tables used by MICM's analytical fixtures."""
import ast
import json
import re
from pathlib import Path

root = Path(__file__).resolve().parents[1]
source = (root.parent/'micm/test/integration/analytical_policy.hpp').read_text()
tables = {}
for name in ('Robertson','Oregonator','Hires','E5'):
    start = source.index('\n{',source.index('void TestAnalytical'+name))
    start = source.index('analytical_concentrations = {',start)
    start = source.index('{',start)
    level = 0
    for end in range(start,len(source)):
        level += (source[end] == '{')-(source[end] == '}')
        if not level:
            break
    literal = re.sub(r'//[^\n]*','',source[start:end+1])
    tables[name] = ast.literal_eval(literal.replace('{','[').replace('}',']'))
(root/'tests/analytical-tables.json').write_text(json.dumps(tables,indent=2)+'\n')
