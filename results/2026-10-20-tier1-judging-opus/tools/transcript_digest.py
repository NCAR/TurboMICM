#!/usr/bin/env python3
"""Condense a Codex session transcript into a timeline (times in MDT, UTC-6).

    transcript_digest.py RUN OUT.txt

Reads transcripts/RUN.jsonl.gz from TurboMICM's `transcripts` branch. Writes one
line per assistant message (commentary and final), per shell command, and per
output line that reports a time, a test count, or a gate result.
"""
import gzip, io, json, re, subprocess, sys
from pathlib import Path
from datetime import datetime, timedelta

run, out = sys.argv[1], sys.argv[2]
blob = subprocess.run(["git", "-C", str(Path.home() / "NSF-NCAR/TurboMICM"), "show", f"transcripts:transcripts/{run}.jsonl.gz"],
                      capture_output=True, check=True).stdout
KEEP = re.compile(r"(elapsed_ms|\bms\b|PASS|FAIL|N/A|passed|failed|max(imum)? rel|largest|speedup|faster|slower|gate|steps=|rejected|median)", re.I)
CMD = re.compile(r'cmd\s*:\s*"((?:[^"\\]|\\.)*)"')

def mdt(ts):
    t = datetime.fromisoformat(ts.replace("Z", "+00:00")) - timedelta(hours=6)
    return t.strftime("%m-%d %H:%M:%S")

lines = []
with gzip.open(io.BytesIO(blob), "rt") as f:
    for raw in f:
        o = json.loads(raw)
        p = o.get("payload") or {}
        ts = o.get("timestamp")
        if not ts:
            continue
        kind = p.get("type")
        if kind == "message" and p.get("role") == "assistant":
            text = " ".join(c.get("text", "") for c in p.get("content", []) if isinstance(c, dict))
            text = " ".join(text.split())
            lines.append(f"[{mdt(ts)}] SAY: {text[:1200]}")
        elif kind in ("custom_tool_call", "function_call"):
            inp = p.get("input") or p.get("arguments") or ""
            cmds = CMD.findall(inp)
            if cmds:
                for c in cmds:
                    c = bytes(c, "utf-8").decode("unicode_escape", errors="replace")
                    lines.append(f"[{mdt(ts)}] RUN: {' '.join(c.split())[:300]}")
            else:
                # apply_patch and similar: keep the first line
                first = " ".join(inp.split())[:200]
                lines.append(f"[{mdt(ts)}] TOOL: {first}")
        elif kind in ("custom_tool_call_output", "function_call_output"):
            outp = p.get("output")
            text = ""
            if isinstance(outp, list):
                text = "\n".join(c.get("text", "") for c in outp if isinstance(c, dict))
            elif isinstance(outp, str):
                text = outp
            text = text.replace("\\n", "\n")
            kept = [l.strip() for l in text.split("\n") if KEEP.search(l) and len(l.strip()) < 400]
            for l in kept[:12]:
                lines.append(f"[{mdt(ts)}]   OUT: {l[:300]}")
        elif o.get("type") == "event_msg" and kind in ("task_complete", "task_started", "thread_goal_updated"):
            lines.append(f"[{mdt(ts)}] EVENT: {kind}")
        elif o.get("type") == "compacted":
            lines.append(f"[{mdt(ts)}] EVENT: context compacted")
open(out, "w").write("\n".join(lines) + "\n")
print(f"{run}: {len(lines)} lines -> {out}")
