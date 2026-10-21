#!/usr/bin/env python3
"""Function extents in C by brace depth, after removing comments, strings and preprocessor lines."""
import re

def strip(text):
    out = []
    i, n = 0, len(text)
    while i < n:
        c = text[i]
        if text.startswith("/*", i):
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            out.append("".join("\n" if ch == "\n" else " " for ch in text[i:j]))
            i = j
        elif text.startswith("//", i):
            j = text.find("\n", i)
            j = n if j < 0 else j
            out.append(" " * (j - i))
            i = j
        elif c in "\"'":
            j = i + 1
            while j < n and text[j] != c:
                j += 2 if text[j] == "\\" else 1
            out.append(c + " " * (min(j, n) - i - 1) + c)
            i = j + 1
        else:
            out.append(c)
            i += 1
    s = "".join(out)
    # blank preprocessor lines (with continuations)
    lines = s.split("\n")
    cont = False
    for k, l in enumerate(lines):
        if cont or l.lstrip().startswith("#"):
            cont = l.rstrip().endswith("\\")
            lines[k] = ""
    return "\n".join(lines)

SIG = re.compile(r"([A-Za-z_]\w*)\s*\((?:[^()]|\([^()]*\))*\)\s*$", re.S)

def functions(text):
    s = strip(text)
    res = []
    depth = 0
    pending = []
    line = 1
    start = None
    name = None
    for ch in s:
        if ch == "\n":
            line += 1
        if depth == 0:
            if ch == "{":
                sig = "".join(pending)
                m = SIG.search(sig)
                if m and "=" not in sig and m.group(1) not in ("if", "for", "while", "switch"):
                    name = m.group(1)
                    start = line - sig[sig.rfind(m.group(1)):].count("\n")
                else:
                    name = None
                depth = 1
                pending = []
                continue
            if ch in ";}":
                pending = []
            else:
                pending.append(ch)
        else:
            if ch == "{":
                depth += 1
            elif ch == "}":
                depth -= 1
                if depth == 0:
                    if name:
                        res.append((name, start, line - start + 1))
                    name = None
                    pending = []
    return res

def comment_line_share(text):
    """Share of non-blank lines that carry a comment (strings left alone)."""
    i, n = 0, len(text)
    marks = set()
    line = 0
    while i < n:
        c = text[i]
        if c == "\n":
            line += 1
            i += 1
        elif text.startswith("/*", i):
            j = text.find("*/", i + 2)
            j = n if j < 0 else j + 2
            for k in range(text.count("\n", i, j) + 1):
                marks.add(line + k)
            line += text.count("\n", i, j)
            i = j
        elif text.startswith("//", i):
            marks.add(line)
            j = text.find("\n", i)
            i = n if j < 0 else j
        elif c in "\"'":
            j = i + 1
            while j < n and text[j] != c and text[j] != "\n":
                j += 2 if text[j] == "\\" else 1
            i = j + 1
        else:
            i += 1
    raw = text.split("\n")
    nonblank = [k for k, l in enumerate(raw) if l.strip()]
    return sum(1 for k in nonblank if k in marks), len(nonblank)
