"""Independent analytic linear-system order, conservation, and grouping checks."""
import csv
import json
import math
import subprocess
from pathlib import Path

import numpy as np
from scipy.linalg import expm

root = Path(__file__).resolve().parents[1]
exact = expm(3*np.array([[-2, 0, .1], [2, -.7, 0], [0, .7, -.1]])) @ np.array([.8, .15, .05])
orders = [2, 3, 4, 3, 4, 2, 3, 4, 2, 3]
rows = []
for method, order in enumerate(orders):
    errors = []
    for h in [.2, .1, .05, .025]:
        command = [str(root/'bin/adaptive_probe'), str(method), '1', '0', '0', '1', str(h)]
        c = json.loads(subprocess.check_output(command, text=True))
        y = np.array(c['y']); error = max(abs(y-exact)); errors.append(error)
        assert abs(sum(y)-1) < 2e-13
        assert c['counts'][2] == 0, c
        grouped = json.loads(subprocess.check_output(command+['7'], text=True))
        np.testing.assert_allclose(np.array(grouped['y']).reshape(3, 7),
                                   np.tile(y[:, None], (1, 7)), rtol=5e-14, atol=5e-16)
        assert grouped['counts'] == c['counts'], (method, grouped['counts'], c['counts'])
        rows.append(dict(method=method, h=h, expected_order=order, error=error,
                         measured_order='' if len(errors)<2 else math.log2(errors[-2]/error)))
    observed = math.log2(errors[-2]/errors[-1])
    assert observed > order-.25, (method, order, observed)
with (root/'results/convergence.csv').open('w', newline='') as f:
    writer=csv.DictWriter(f, fieldnames=rows[0]); writer.writeheader(); writer.writerows(rows)
print('10 methods: analytic convergence orders, conservation, and 7-cell grouping pass')
