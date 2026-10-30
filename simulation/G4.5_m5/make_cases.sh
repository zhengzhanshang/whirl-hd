#!/usr/bin/env bash
# =====================================================================
# Build the seven G=4.5, m=5 case folders from input_template.json and
# pairs.txt, so each folder holds what run.sh expects:
#   input.json, run.sh, sw_model.py, spherical-2d-hd-s2s4.py, evp_spherical_hd.py
#
# Added for this repository. It was not used on the cluster; the original
# case folders were set up by hand there. It only fills in s2 and s4, which
# are the only keys in which the seven input.json files differ.
#
#   bash make_cases.sh            # creates ./s2=..._s4=.../ next to this script
#   OUT=/scratch/me bash make_cases.sh
# =====================================================================
set -euo pipefail
HERE="$(cd "$(dirname "$0")" && pwd)"
OUT=${OUT:-$HERE}

grep -v '^#' "$HERE/pairs.txt" | while read -r S2 S4 NAME; do
  [ -n "${NAME:-}" ] || continue
  D="$OUT/$NAME"
  if [ -e "$D" ]; then echo "exists, skipped: $D"; continue; fi
  mkdir -p "$D"
  cp "$HERE"/{run.sh,sw_model.py,spherical-2d-hd-s2s4.py,evp_spherical_hd.py} "$D"/
  python3 - "$HERE/input_template.json" "$D/input.json" "$S2" "$S4" <<'PY'
import json, sys
p = json.load(open(sys.argv[1]))
p["s2"], p["s4"] = float(sys.argv[3]), float(sys.argv[4])
json.dump(p, open(sys.argv[2], "w"), indent=2)
PY
  echo "created $D  (s2=$S2, s4=$S4)"
done
