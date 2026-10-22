#!/bin/zsh
# Close, collect, and archive a finished Tier 2 run.
#
#   admin/collect_tier2.sh RUN ROLLOUT.jsonl
#
# For a run that admin/run_state.py reports DONE:
# 1. quits Codex in its tmux session (socket turbomicm, session tier2-RUN) with
#    /quit, and closes the session;
# 2. checks that the run's working tree is clean (ignored files aside);
# 3. fetches its branch main into this repository's tier2/<model>/<effort> and
#    pushes it;
# 4. archives the session record with admin/archive_transcript.py as
#    transcripts/tier2/RUN.jsonl.gz on the transcripts branch, from a temporary
#    worktree, adds its row to the branch README's Tier 2 table, scans it for
#    credential patterns (reporting the fields the matches are in), commits,
#    and pushes.
set -eu
here=${0:A:h}
repo=${here:h}
run=$1
record=${2:A}
model=${run%%-*}
effort=${run#*-}
dir=$HOME/NSF-NCAR/ODE/benchmarks/tier2/$run
tmux_bin=""
for t in /opt/homebrew/bin/tmux $HOME/homebrew/bin/tmux; do [[ -x $t ]] && { tmux_bin=$t; break; }; done
session=${${record:t:r}[-36,-1]}  # rollout-<time>-<session id>.jsonl

if $tmux_bin -L turbomicm has-session -t tier2-$run 2>/dev/null; then
  $tmux_bin -L turbomicm send-keys -t tier2-$run -l "/quit"
  sleep 1
  $tmux_bin -L turbomicm send-keys -t tier2-$run Enter
  sleep 5
  echo "quit: $(date '+%H:%M:%S')"
  $tmux_bin -L turbomicm kill-session -t tier2-$run
fi
dirty=$(git -C $dir --no-optional-locks status --porcelain)
[[ -z $dirty ]] || { echo "STOP: $run's working tree is not clean:"; echo $dirty; exit 1; }
git -C $repo fetch -q $dir main:tier2/$model/$effort
echo "collected: tier2/$model/$effort at $(git -C $repo log --oneline -1 tier2/$model/$effort)"
git -C $repo push -q origin tier2/$model/$effort 2>&1 | command grep -v '^remote:' || true

work=$(mktemp -d)
git -C $repo fetch -q origin
git -C $repo worktree add -q --detach $work/transcripts origin/transcripts
cd $work/transcripts
mkdir -p transcripts/tier2
python3 $repo/admin/archive_transcript.py $record transcripts/tier2/$run.jsonl.gz
gunzip -c transcripts/tier2/$run.jsonl.gz | python3 -c "
import json, re, sys
pattern = re.compile(r'(sk-[A-Za-z0-9]{20}|ghp_[A-Za-z0-9]{20}|gho_[A-Za-z0-9]{20}|github_pat_[A-Za-z0-9_]{20})')
hits = {}
for line in sys.stdin:
    item = json.loads(line)
    payload = item.get('payload') if isinstance(item.get('payload'), dict) else {}
    for key, value in payload.items():
        text = value if isinstance(value, str) else json.dumps(value)
        for _ in pattern.findall(text):
            hits[key] = hits.get(key, 0) + 1
print('credential-pattern matches by field:', hits or 'none')
"
size=$(python3 -c "import os, sys; print(f'{os.path.getsize(sys.argv[1]) / 1e6:.1f} MB')" $record)
python3 - $run $session "$size" <<'EOF'
import sys
from pathlib import Path
run, session, size = sys.argv[1:]
readme = Path("README.md")
lines = readme.read_text().splitlines(keepends=True)
rows = [i for i, line in enumerate(lines) if line.startswith("| `transcripts/tier2/")]
row = f"| `transcripts/tier2/{run}.jsonl.gz` | Tier 2 {run} | Codex CLI 0.163.0 | `{session}` | {size} |\n"
if any(f"transcripts/tier2/{run}.jsonl.gz" in line for line in lines):
    sys.exit(f"README already lists {run}")
lines.insert(rows[-1] + 1, row)
readme.write_text("".join(lines))
print("README:", row.strip())
EOF
git add README.md transcripts/tier2/$run.jsonl.gz
git commit -q -m "Tier 2 $run transcript (session ${session[1,8]})

Co-Authored-By: Claude Opus 5.5 <noreply@anthropic.com>"
git push -q origin HEAD:transcripts
echo "archived: $(git log --oneline -1)"
cd $repo
git worktree remove $work/transcripts
rm -rf $work
