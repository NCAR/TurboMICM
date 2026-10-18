#!/bin/zsh
# Launch a Codex run as the SSE launches one by hand, from a script.
#
#   admin/launch_codex_run.sh [--tier3] RUN MODEL EFFORT [--check-only]
#
# For example: admin/launch_codex_run.sh luna-high gpt-6-luna high
#              admin/launch_codex_run.sh --tier3 sol-xhigh gpt-6.1-sol xhigh
#
# 1. Checks, stopping without launching if any fails:
#    - Tier 1: the run directory is a git repository with no commits whose only
#      file is PROMPT.md, identical to prompts/RUN.txt;
#    - Tier 3: the run directory is an untouched clone of its Tier 1 run, made by
#      tier3/admin/setup_runs.py: HEAD is the Tier 1 commit the prompt names, no
#      remote, and nothing changed or added but PROMPT.md, identical to
#      prompts/RUN.txt;
#    - no CLAUDE.md or AGENTS.md in the run root, above it, or in ~/.codex; no
#      other Codex session running.
# 2. Starts a detached tmux server (socket "turbomicm", one session named after the
#    run, tier3-RUN for Tier 3) with an empty environment, so nothing from the
#    calling program reaches Codex; the session runs the SSE's login shell, which
#    sets up the SSE's usual environment, in the run directory.
# 3. Types the standard command (Tier 3 adds --search, for the literature), answers
#    the folder-trust prompt, waits for the composer, types the standard goal,
#    checks it on screen, and submits it.
# 4. Finds the new session record and confirms its working directory, Codex
#    version, base instructions, model, effort, and goal, and the Codex process's
#    command line and environment.
# On any screen it does not expect, it stops and prints the screen. Watch or take
# over the session with: tmux -L turbomicm attach -t RUN (tier3-RUN for Tier 3)
#
# LAUNCH_ROOT replaces the run root, for a test on a throwaway root.
set -u
tier=1
[[ ${1:-} == --tier3 ]] && { tier=3; shift; }
run=$1 model=$2 effort=$3 mode=${4:-launch}
if (( tier == 3 )); then
  root=${LAUNCH_ROOT:-$HOME/NSF-NCAR/ODE/benchmarks/tier3}
  session=tier3-$run
  command="codex --search -m $model -c model_reasoning_effort=$effort -c features.memories=false"
else
  root=${LAUNCH_ROOT:-$HOME/NSF-NCAR/ODE/benchmarks}
  session=$run
  command="codex -m $model -c model_reasoning_effort=$effort -c features.memories=false"
fi
dir=$root/$run
goal='/goal read PROMPT.md and work until done, commit'
tmux_bin=/opt/homebrew/bin/tmux
tm() { $tmux_bin -L turbomicm "$@"; }
screen() { tm capture-pane -p -t $session 2>/dev/null | sed -e '/^[[:space:]]*$/d'; }
stop() { echo "STOP: $1"; [[ ${2:-} == show ]] && screen | tail -20; exit 1; }
wait_for() {  # wait_for TEXT SECONDS: true once TEXT is on screen
  local i=0
  while (( i < $2 )); do screen | command grep -q -F -- "$1" && return 0; sleep 1; (( i++ )); done
  return 1
}
codex_tuis() { ps -axo pid=,command= | awk '$2 == "codex" {print $1}'; }

echo "check: $(codex --version)"
[[ -d $dir/.git ]] || stop "$dir is not a git repository"
if (( tier == 3 )); then
  tier1=$(tr '\n' ' ' < $root/prompts/$run.txt | command grep -o -E 'repository +at +commit +[0-9a-f]{7,40}' | awk '{print $NF}')
  [[ -n $tier1 ]] || stop "prompts/$run.txt names no Tier 1 commit"
  [[ "$(git -C $dir rev-parse HEAD)" == "$(git -C $dir rev-parse -q --verify "$tier1^{commit}")" ]] \
    || stop "$run is not at its Tier 1 commit $tier1"
  [[ -z "$(git -C $dir remote)" ]] || stop "$run has a remote"
  # PROMPT.md is tracked in some Tier 1 repositories (" M") and ignored in others ("!!")
  extra=$(git --no-optional-locks -C $dir status --porcelain --ignored=matching | command grep -v -x -E '( M|!!) PROMPT.md')
  [[ -z $extra ]] || stop "$run has changes besides PROMPT.md: $(echo $extra | tr '\n' ' ')"
  state="untouched at Tier 1 commit $tier1"
else
  git -C $dir rev-parse -q --verify HEAD >/dev/null && stop "$run already has commits"
  [[ "$(command ls -A $dir | command grep -v -x '.git')" == PROMPT.md ]] || stop "$run holds more than PROMPT.md"
  state="untouched"
fi
cmp -s $dir/PROMPT.md $root/prompts/$run.txt || stop "PROMPT.md differs from prompts/$run.txt"
d=$root
while [[ $d != / ]]; do
  [[ -e $d/CLAUDE.md || -e $d/AGENTS.md ]] && stop "an instruction file is in $d"
  d=$d:h
done
[[ -e $HOME/.codex/CLAUDE.md || -e $HOME/.codex/AGENTS.md ]] && stop "an instruction file is in $HOME/.codex"
[[ -n "$(codex_tuis)" ]] && stop "a Codex session is running: $(codex_tuis | tr '\n' ' ')"
tm has-session -t $session 2>/dev/null && stop "tmux session $session already exists"
echo "check: $run is $state, PROMPT.md matches prompts/$run.txt, no instruction files, no Codex session running"
[[ $mode == --check-only ]] && exit 0

marker=$(mktemp)
env -i HOME=$HOME USER=$USER LOGNAME=$USER SHELL=/bin/zsh TERM=xterm-256color LANG=en_US.UTF-8 \
    PATH=/usr/bin:/bin:/usr/sbin:/sbin $tmux_bin -L turbomicm new-session -d -s $session -x 200 -y 50 -c $dir \
  || stop "tmux did not start"
tm send-keys -t $session -l "$command"
tm send-keys -t $session Enter
echo "launch: $(date '+%Y-%m-%d %H:%M:%S') $command"

if wait_for 'Trust this folder?' 30; then
  screen | command grep -q -F '› 1. Trust and continue' || stop "unexpected trust prompt" show
  tm send-keys -t $session Enter
  echo "launch: trusted the run directory"
fi
wait_for 'Ask Codex to do anything' 30 || stop "the composer did not appear" show
screen | command grep -q -i -F -- "$model" || stop "$model is not on the status line" show
tm send-keys -t $session -l "$goal"
wait_for "› $goal" 10 || stop "the goal did not appear in the composer" show
tm send-keys -t $session Enter
echo "launch: $(date '+%Y-%m-%d %H:%M:%S') sent the goal"

record=""
for i in {1..60}; do
  for f in $(find $HOME/.codex/sessions -name 'rollout-*.jsonl' -newer $marker 2>/dev/null); do
    head -c 20000 $f | command grep -q -F "\"cwd\":\"$dir\"" && command grep -q -F '"thread_goal_updated"' $f && record=$f
  done
  [[ -n $record ]] && break
  sleep 1
done
rm -f $marker
[[ -n $record ]] || stop "no session record with the goal after 60 s" show
pid=""
for p in $(codex_tuis); do
  [[ "$(lsof -a -p $p -d cwd -Fn 2>/dev/null | command grep '^n' | cut -c2-)" == $dir ]] && pid=$p
done
echo "record: $record"
echo "pid: $pid"
python3 - "$record" "$dir" "$model" "$effort" <<'EOF'
import json, sys
path, cwd, model, effort = sys.argv[1:]
lines = [json.loads(l) for l in open(path)]
meta = lines[0]['payload']
settings, goal, switch = set(), None, False
for d in lines:
    p = d.get('payload') if isinstance(d.get('payload'), dict) else {}
    if p.get('type') == 'thread_settings_applied':
        settings.add((p['thread_settings'].get('model'), p['thread_settings'].get('reasoning_effort')))
    if d.get('type') == 'turn_context':
        settings.add((p.get('model'), p.get('effort')))
    if p.get('type') == 'thread_goal_updated' and goal is None:
        goal = p['goal'].get('objective')
    switch = switch or '<model_switch>' in json.dumps(p)
base = meta['base_instructions']['provenance'].get('model')
ok = meta['cwd'] == cwd and base == model and settings == {(model, effort)} and not switch \
     and goal == 'read PROMPT.md and work until done, commit'
print(f"verify: cwd {meta['cwd']}; Codex {meta['cli_version']}; base instructions {base}; "
      f"model and effort {sorted(settings)}; model switch {switch}; goal {goal!r}: {'OK' if ok else 'MISMATCH'}")
sys.exit(0 if ok else 1)
EOF
verified=$?
if [[ -n $pid ]]; then
  echo "verify: command: $(ps -o command= -p $pid)"
  echo "verify: Claude variables in its environment: $(ps -wwE -p $pid -o command= | tr ' ' '\n' | command grep -E '^(CLAUDE|AI_AGENT|ANTHROPIC)' | cut -d= -f1 | tr '\n' ' ')"
fi
exit $verified
