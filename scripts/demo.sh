#!/usr/bin/env bash
# End-to-end demo against a running service (default http://localhost:8000).
# Shows: upload -> progress polling -> download verification, amplicon merge,
# plus cancellation and a rejected (malformed) upload.
set -euo pipefail

BASE_URL="${BASE_URL:-http://localhost:8000}"
WORK="$(mktemp -d)"
trap 'rm -rf "$WORK"' EXIT

need() { command -v "$1" >/dev/null || { echo "missing tool: $1" >&2; exit 1; }; }
need curl
need python3

echo "== health =="
curl -fsS "$BASE_URL/health"
echo

mk() { # file id seq side
  printf '@%s/%s desc\n%s\n+\n%s\n' "$2" "$4" "$3" "$(printf '~%.0s' $(seq ${#3}))" >> "$1"
}

echo "== generating shuffled R1/R2 =="
R1="$WORK/R1.fastq"; R2="$WORK/R2.fastq"
for id in alpha beta gamma delta epsilon; do
  mk "$R1" "$id" "ACGTACGTAC" 1
done
for id in epsilon gamma alpha delta beta; do   # deliberately shuffled
  mk "$R2" "$id" "TGCATGCATG" 2
done

echo "== upload =="
TASK=$(curl -fsS -X POST "$BASE_URL/tasks" \
  -F "r1=@$R1;filename=R1.fastq" \
  -F "r2=@$R2;filename=R2.fastq")
echo "$TASK"
ID=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$TASK")

echo "== poll until terminal =="
while :; do
  BODY=$(curl -fsS "$BASE_URL/tasks/$ID")
  STATUS=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["status"])' <<<"$BODY")
  echo "$BODY"
  case "$STATUS" in
    succeeded|failed|cancelled) break ;;
  esac
  sleep 0.5
done

echo "== download and verify order =="
curl -fsS "$BASE_URL/tasks/$ID/download/r1" -o "$WORK/out_R1.fastq"
curl -fsS "$BASE_URL/tasks/$ID/download/r2" -o "$WORK/out_R2.fastq"
python3 - "$WORK/out_R1.fastq" "$WORK/out_R2.fastq" <<'PY'
import sys
ids = []
for path in sys.argv[1:]:
    got = []
    with open(path, "rb") as fh:
        while True:
            header = fh.readline()
            if not header:
                break
            got.append(header[1:].split()[0].decode().rsplit("/", 1)[0])
            fh.readline(); fh.readline(); fh.readline()
    ids.append(got)
print("R1 ids:", ids[0])
print("R2 ids:", ids[1])
assert ids[0] == ids[1] == ["alpha", "beta", "gamma", "delta", "epsilon"]
print("OK: paired in R1 order")
PY

echo "== amplicon merge on the succeeded task =="
# Build a separate task whose 60bp mates share a 30bp overlap so the merge
# endpoint has an admissible, unique-longest candidate to stitch.
MR1="$WORK/merge_R1.fastq"; MR2="$WORK/merge_R2.fastq"
python3 - "$MR1" "$MR2" <<'PY'
import sys
comp = bytes.maketrans(b"ACGT", b"TGCA")
def rec(rid, seq, qual, suf):
    return b"@" + rid.encode() + suf + b"\n" + seq + b"\n+\n" + qual + b"\n"
# rc(R2) = C*30 + G*30 aligns after the 30 A bases of R1 -> 30bp overlap.
r1 = rec("merge-1", b"A" * 30 + b"C" * 30, b"E" * 60, b"/1")
r2rc = b"C" * 30 + b"G" * 30
r2 = rec("merge-1", r2rc.translate(comp)[::-1], b"F" * 60, b"/2")
# A second pair with unrelated random reads -> no candidate -> skipped.
r1 += rec("merge-2", b"ACGTACGGTA" * 4, b"H" * 40, b"")
r2 += rec("merge-2", b"TGCATGCATA" * 4, b"I" * 40, b"")
open(sys.argv[1], "wb").write(r1)
open(sys.argv[2], "wb").write(r2)
PY
MTASK=$(curl -fsS -X POST "$BASE_URL/tasks" \
  -F "r1=@$MR1;filename=merge_R1.fastq" \
  -F "r2=@$MR2;filename=merge_R2.fastq")
MID=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$MTASK")
while :; do
  ST=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["status"])' \
    < <(curl -fsS "$BASE_URL/tasks/$MID"))
  [ "$ST" = "succeeded" ] && break
  sleep 0.3
done
curl -fsS -X POST "$BASE_URL/tasks/$MID/merge" | python3 -m json.tool

echo "== malformed input is rejected as a whole task =="
printf '@x\nACGT\n+\n~~\n' > "$WORK/bad_R2.fastq"     # seq 4 != qual 2
BAD=$(curl -fsS -X POST "$BASE_URL/tasks" \
  -F "r1=@$R1;filename=R1.fastq" \
  -F "r2=@$WORK/bad_R2.fastq;filename=R2.fastq")
BAD_ID=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$BAD")
sleep 1
curl -fsS "$BASE_URL/tasks/$BAD_ID"
echo
echo "(download attempt on failed task -> HTTP 409)"
curl -s -o /dev/null -w '%{http_code}\n' "$BASE_URL/tasks/$BAD_ID/download/r1"

echo "== cancellation returns immediately, task never publishes =="
BIG_R1="$WORK/big_R1.fastq"; BIG_R2="$WORK/big_R2.fastq"
python3 - "$BIG_R1" "$BIG_R2" <<'PY'
import sys
n = 200_000
with open(sys.argv[1], "wb") as f1, open(sys.argv[2], "wb") as f2:
    for i in range(n):
        f1.write(b"@r%09d/1\n" % i + b"A" * 40 + b"\n+\n" + b"~" * 40 + b"\n")
    for i in range(n - 1, -1, -1):
        f2.write(b"@r%09d/2\n" % i + b"T" * 40 + b"\n+\n" + b"!" * 40 + b"\n")
PY
CANCEL=$(curl -fsS -X POST "$BASE_URL/tasks" \
  -F "r1=@$BIG_R1;filename=big_R1.fastq" \
  -F "r2=@$BIG_R2;filename=big_R2.fastq")
CANCEL_ID=$(python3 -c 'import json,sys; print(json.load(sys.stdin)["id"])' <<<"$CANCEL")
curl -fsS -X POST "$BASE_URL/tasks/$CANCEL_ID/cancel"
echo
sleep 2
curl -fsS "$BASE_URL/tasks/$CANCEL_ID"
echo

echo "== cancelling the already succeeded task -> HTTP 409 =="
curl -s -o /dev/null -w '%{http_code}\n' -X POST "$BASE_URL/tasks/$ID/cancel"

echo "ALL DEMO STEPS DONE"
