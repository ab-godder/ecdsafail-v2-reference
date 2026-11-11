#!/bin/bash
# Reproduces every receipt in receipts/ outside baseline/, and the tracked artifacts/.
# The harness steps run only when ECDSAFAIL_V2 names the benchmark checkout, which
# harness/Cargo.toml reads (read-only) from ../ecdsafail-v2 next to this repo.
set -euo pipefail
cd "$(dirname "$0")/.."
PY=${PY:-python3}
export PATH="$HOME/.cargo/bin:$PATH"
exec > >(tee receipts/run_all.log) 2>&1

# BSD time takes -l and GNU time takes -v. Without either, the steps run untimed.
TIME_FLAG=
for flag in -l -v; do
  if /usr/bin/time "$flag" true > /dev/null 2>&1; then TIME_FLAG=$flag; break; fi
done
timed() {
  if [ -n "$TIME_FLAG" ]; then /usr/bin/time "$TIME_FLAG" "$@"; else "$@"; fi
}
timed_to() {
  local out=$1; shift
  if [ -n "$TIME_FLAG" ]; then /usr/bin/time "$TIME_FLAG" -o "$out" "$@"; else echo "untimed: no /usr/bin/time with -l or -v" > "$out"; "$@"; fi
}

$PY -I scripts/model_check.py > receipts/model_check.jsonl
$PY -I scripts/width_check.py > receipts/width_check.jsonl
$PY -I scripts/check.py --toy-bits 4 6 7 8 9 10 --toy-per-size 6 --out receipts/toy_exhaustive.json > /dev/null
$PY -I scripts/mutation_control.py > receipts/mutation_control.json
$PY -I scripts/component_check.py > receipts/secp256k1_components.json
timed $PY -I scripts/check.py --secp --secp-random 120 --seed 2026 --out receipts/secp256k1_sampled.json > receipts/secp256k1_sampled.log 2>&1
$PY -I scripts/breakdown.py > receipts/secp256k1_breakdown.json

$PY -I scripts/emit.py --curve toy:127 --addend-index 5 --out artifacts/toy127_circuit.json
$PY -I scripts/emit.py --out artifacts/secp256k1_G_circuit.json.gz
$PY -I scripts/emit.py --curve toy:127 --addend-index 5 --hier --out artifacts/toy127_circuit.hier.json
$PY -I scripts/emit.py --hier --out artifacts/secp256k1_G_circuit.hier.json
$PY -I scripts/make_cases.py --curve toy:127 --addend-index 5 --out artifacts/toy127_cases.json
$PY -I scripts/make_cases.py --random 250 --out artifacts/secp256k1_cases.json

if [ -z "${ECDSAFAIL_V2:-}" ]; then
  echo "skipping the harness steps: set ECDSAFAIL_V2 to the ecdsafail-v2 benchmark checkout to run them"
  echo "receipts/harness_*.json and artifacts/toy127_submission.json keep their committed contents"
elif [ ! -f "$ECDSAFAIL_V2/contracts/ir/Cargo.toml" ]; then
  echo "ECDSAFAIL_V2=$ECDSAFAIL_V2 is not an ecdsafail-v2 checkout (no contracts/ir/Cargo.toml)" >&2
  exit 1
elif [ "$(cd "$ECDSAFAIL_V2" && pwd -P)" != "$(cd ../ecdsafail-v2 2>/dev/null && pwd -P)" ]; then
  echo "harness/Cargo.toml reads the benchmark from $(cd .. && pwd)/ecdsafail-v2, which is not $ECDSAFAIL_V2" >&2
  echo "clone or symlink the benchmark checkout there" >&2
  exit 1
else
  (cd harness && CARGO_TARGET_DIR="$PWD/target" cargo build --release --offline --quiet)
  H=harness/target/release/qec-harness
  $H --circuit artifacts/toy127_circuit.json --cases artifacts/toy127_cases.json --lower --write-submission artifacts/toy127_submission.json > receipts/harness_toy127.json
  gzip -dc artifacts/secp256k1_G_circuit.json.gz | timed_to receipts/harness_secp256k1.time $H --circuit - --cases artifacts/secp256k1_cases.json --lower > receipts/harness_secp256k1.json
fi
echo "all receipts written in ${SECONDS}s"
