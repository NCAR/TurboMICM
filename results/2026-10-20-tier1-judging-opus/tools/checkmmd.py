#!/usr/bin/env python3
"""Render every ```mermaid block of a Markdown file with mmdc; report failures.

    checkmmd.py PAGE.md [PAGE.md ...]
"""
import re, subprocess, sys, tempfile
from pathlib import Path

bad = 0
for page in sys.argv[1:]:
    text = Path(page).read_text()
    blocks = re.findall(r"```mermaid\n(.*?)```", text, re.S)
    for i, b in enumerate(blocks, 1):
        with tempfile.TemporaryDirectory() as d:
            src, out = Path(d) / "d.mmd", Path(d) / "d.svg"
            src.write_text(b)
            r = subprocess.run(["mmdc", "-q", "-i", str(src), "-o", str(out)], capture_output=True, text=True)
            ok = r.returncode == 0 and out.exists() and "Syntax error" not in out.read_text(errors="replace")
            first = b.strip().split("\n")[0]
            print(f"{page}: diagram {i} ({first}): {'ok' if ok else 'FAILED'}")
            if not ok:
                bad += 1
                print((r.stderr or r.stdout)[-1500:])
sys.exit(1 if bad else 0)
