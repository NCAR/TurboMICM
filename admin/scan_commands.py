#!/usr/bin/env python3
"""List a Codex run's commands that reach outside its directory, install, or download.

    admin/scan_commands.py ROLLOUT.jsonl RUN_DIR [--allow PATH ...]

Reads the session record (read-only) and prints, for each kind of finding, how
many commands match and the first few, with their times:
- paths outside RUN_DIR: absolute paths, ~ paths, and ../ paths that leave it,
  other than the allowed inputs (by default RUN_DIR/../micm, ../reference, and
  ../scenarios, Codex's own AGENTS.md checks in the parent folders, the system's
  and Homebrew's tools, conda, and /tmp);
- other runs and the experiment's private folders: any mention of a sibling
  run directory (<model>-<effort>), TurboMICM, tier2-work, .tier2-heldout, or
  the transcripts;
- installs and downloads: brew, pip, conda or mamba install, curl, wget, git
  clone or fetch, npm, and web access;
- connector, MCP, and web-search calls, which are not shell commands: Codex
  0.169.0 records them as McpToolCall items, function calls in an mcp__
  namespace, or code-mode JavaScript naming an mcp__ tool (ChatGPT's app
  connectors are mcp__codex_apps__*);
- subagents: spawn and message calls, and messages between agents (each
  subagent writes a record of its own, to be scanned and counted too).
A finding is a lead for review, not a verdict: the command's context decides.
Needs only the Python standard library.
"""

import argparse
import json
import re
from pathlib import Path

parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
parser.add_argument("record", type=Path)
parser.add_argument("run_dir", type=Path)
parser.add_argument("--allow", nargs="*", default=[], help="more allowed path prefixes")
parser.add_argument("--show", type=int, default=6, help="examples per finding (default 6)")
args = parser.parse_args()

run = str(args.run_dir.resolve())
root = str(args.run_dir.resolve().parent)
home = str(Path.home())
inputs = [Path(root) / name for name in ("micm", "reference", "scenarios")]
allowed = [run, *map(str, inputs), *(str(p.resolve()) for p in inputs), "/usr/", "/bin/", "/sbin/", "/opt/",
           "/System/", "/Library/Developer/", "/private/tmp", "/tmp", "/dev/", f"{home}/homebrew", f"{home}/miniforge3",
           f"{home}/.codex/", f"{home}/Library/Caches", *args.allow]
agents_check = re.compile(r"(AGENTS|CLAUDE)\.md")
siblings = re.compile(r"\b(sonnet|opus|fable|luna|sol|astra|terra)-(high|xhigh|max)\b")
private = re.compile(r"TurboMICM|tier2-work|tier1-work|tier3-work|\.tier[0-9]-heldout|transcripts/|aborted-runs")
install = re.compile(r"\bbrew\s+(install|upgrade|reinstall)|\bpip3?\s+install|\b(conda|mamba)\s+install|\bcurl\b|"
                     r"\bwget\b|\bgit\s+(clone|fetch|pull|submodule)|\bnpm\b|https?://")
path_token = re.compile(r"(?<![\w.$-])(~/[^\s'\";|&)`,\\]*|/[A-Za-z][^\s'\";|&)`,\\]*|(?:\.\./)+[^\s'\";|&)`,\\]*)")


def commands(record):
    """(time, command text) for every shell command in the record."""
    with open(record, errors="replace") as file:
        for line in file:
            try:
                item = json.loads(line)
            except ValueError:
                continue
            payload = item.get("payload") if isinstance(item.get("payload"), dict) else {}
            if payload.get("type") == "item_completed":
                command = (payload.get("item") or {}).get("command")
                if isinstance(command, list) and command:
                    cwd = ((payload.get("item") or {}).get("cwd") or run).removeprefix("file://")
                    yield item.get("timestamp", "")[11:19], command[-1], cwd


def outside(text, cwd):
    found = []
    for token in path_token.findall(text):
        if token.startswith("~/"):
            path = home + token[1:]
        elif token.startswith("../"):
            path = str((Path(cwd) / token).resolve())
        else:
            path = token
            # an absolute path must start in a folder that exists: "/len(x" or "/name" is code, not a path
            if not Path("/" + path.split("/")[1]).is_dir():
                continue
        if not any(path == a.rstrip("/") or path.startswith(a.rstrip("/") + "/") for a in allowed):
            found.append(token)
    return found


findings = {"paths outside the run directory and its inputs": [], "other runs and private folders": [],
            "installs and downloads": []}
total = 0
for stamp, text, cwd in commands(args.record):
    total += 1
    # Follow "cd DIR" within the command, so that relative paths after it resolve from DIR
    away, here = [], cwd
    for segment in re.split(r"&&|\|\||;|\n", text):
        moved = re.match(r"\s*(?:do\s+|then\s+|else\s+)?\(?\s*cd\s+([^\s;&|()]+)\s*$", segment)
        if moved:
            target = moved.group(1).strip("'\"")
            here = str((Path(here) / target).resolve()) if not target.startswith(("/", "~")) else \
                str(Path(target.replace("~", home, 1)).resolve())
            continue
        # Codex's own check for instruction files in the parent folders names them, and nothing else;
        # a search for them by name (rg --files -g AGENTS.md, find -name AGENTS.md) reads no file
        if re.search(r"\b(rg|find|ls|test)\b.*(AGENTS|CLAUDE)\.md", segment):
            continue
        away += [t for t in outside(segment, here) if not re.search(r"(AGENTS|CLAUDE)\.md$", t)]
    if away:
        findings["paths outside the run directory and its inputs"].append((stamp, text, away))
    if siblings.search(text.replace(Path(run).name, "")) or private.search(text):
        findings["other runs and private folders"].append((stamp, text, []))
    if install.search(text):
        findings["installs and downloads"].append((stamp, text, []))

# Tool calls that are not shell commands (Codex 0.169.0's record): calls to MCP servers and
# ChatGPT's app connectors (an McpToolCall item, a function call in an mcp__ namespace, or
# code-mode JavaScript that names an mcp__ tool), web searches, and messages to subagents
connector = re.compile(r"mcp__(codex_apps|[\w-]+)__\w+|\bweb_search\b")
findings["connector, MCP, and web-search calls"] = []
findings["subagents"] = []
with open(args.record, errors="replace") as file:
    for line in file:
        try:
            item = json.loads(line)
        except ValueError:
            continue
        payload = item.get("payload") if isinstance(item.get("payload"), dict) else {}
        stamp, kind = item.get("timestamp", "")[11:19], payload.get("type")
        inner = payload.get("item") or {}
        if kind == "item_completed" and inner.get("type") in ("McpToolCall", "WebSearch"):
            findings["connector, MCP, and web-search calls"].append(
                (stamp, f"{inner.get('type')} {inner.get('server', '')} {inner.get('tool', '')}", []))
        elif kind == "function_call" and str(payload.get("namespace", "")).startswith("mcp__"):
            findings["connector, MCP, and web-search calls"].append(
                (stamp, f"{payload['namespace']} {payload.get('name')}", []))
        elif kind == "custom_tool_call" and connector.search(payload.get("input", "")):
            names = sorted(set(m.group(0) for m in connector.finditer(payload["input"])))
            findings["connector, MCP, and web-search calls"].append((stamp, f"exec names {' '.join(names)}", []))
        elif kind == "web_search_call":
            findings["connector, MCP, and web-search calls"].append((stamp, "web_search_call", []))
        elif kind == "agent_message" or (kind == "function_call" and payload.get("namespace") == "collaboration"):
            findings["subagents"].append((stamp, f"{kind} {payload.get('name', '')} {payload.get('recipient', '')}", []))

print(f"{args.record.name}: {total} commands")
for kind, rows in findings.items():
    print(f"{kind}: {len(rows)}")
    for stamp, text, away in rows[:args.show]:
        note = f"  [{', '.join(sorted(set(away)))[:160]}]" if away else ""
        print(f"  {stamp} {' '.join(text.split())[:220]}{note}")
