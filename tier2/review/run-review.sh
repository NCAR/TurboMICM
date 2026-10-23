#!/bin/zsh
# Run one reviewer of Tier 2's blind review, in a fresh Claude Code session.
#
#   tier2/review/run-review.sh MODEL ORDER REVIEW_DIR OUT_DIR
#
# For example: tier2/review/run-review.sh opus "A B C D E F G H I" ~/NSF-NCAR/tier2-work/blind-review/REVIEW ~/NSF-NCAR/tier2-work/blind-review
#
# As Tier 3's review did (results/2026-10-18-tier3-blind-review/run-opus.sh): a
# non-interactive session (claude -p) at effort max, started in REVIEW_DIR, which
# holds REVIEW.md, PROMPT.md, the anonymized entries, scenarios/, reference/, and
# mechanisms/; tools limited to Read, Grep, Glob, and Bash restricted to wc, grep,
# head, tail, and cat; no web, no agents, no edits. The report is the session's
# final message, saved with the session record (OUT_DIR/report-MODEL.json) and
# the times (OUT_DIR/report-MODEL.err); OUT_DIR/report-MODEL.done marks the end.
set -u
model=$1
order=$2
review=${3:A}
out=${4:A}
order_words=(${=order})
listed="${(j:, :)order_words}"
prompt="You are one of the two reviewers in a blind review. Your working directory holds REVIEW.md, PROMPT.md, nine entries ($listed), scenarios/, reference/, and mechanisms/. Read REVIEW.md first and follow it exactly. Review the entries in this order: $listed; then compare them. Read only inside this working directory: do not read anything elsewhere on this machine, do not use the web, and do not try to identify the agents. Do not build or run any code; the only shell commands you may use are grep, wc, head, tail, and cat, besides the Read, Grep, and Glob tools. Take the time to read each entry's generator thoroughly, with its tests, build, and documents, and do the three exercises on paper for each. Your final message must be the complete report in Markdown, in the order REVIEW.md gives, ending with the JSON block it specifies. Nothing but that final message is collected, so put the whole report in it."
cd $review || exit 1
echo "start $(date '+%Y-%m-%d %H:%M:%S') model $model order $order" > $out/report-$model.err
env -u CLAUDECODE claude -p --model $model --effort max --tools "Read,Grep,Glob,Bash" \
    --allowedTools "Bash(wc *) Bash(grep *) Bash(head *) Bash(tail *) Bash(cat *)" \
    --disallowedTools "WebSearch WebFetch Agent Edit Write NotebookEdit Bash(git *) Bash(ls *) Bash(find *) Bash(cd *) Bash(python*) Bash(make *) Bash(curl *)" \
    --output-format json "$prompt" > $out/report-$model.json 2>> $out/report-$model.err
echo "exit $? $(date '+%Y-%m-%d %H:%M:%S')" >> $out/report-$model.err
touch $out/report-$model.done
