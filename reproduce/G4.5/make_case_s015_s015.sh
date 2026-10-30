#!/usr/bin/env bash
# =====================================================================
# Create and (optionally) submit the missing G=4.5, m=5 case
#     s2=0.15_s4=0.15
# by cloning an existing m=5 case directory and changing ONLY s2 and s4.
#
# Run from:  paper-redo_G=4.5/m=3-15/m=5_s2s4_sweep/
#
# DRY RUN BY DEFAULT. It builds the directory and shows you the diff, but
# submits nothing. When the diff looks right:
#     SUBMIT=1 bash make_case_s015_s015.sh
#
# Template: s2=0.15_s4=0.17 (closest pair; same s2). Override with
#     TEMPLATE=<dir> bash make_case_s015_s015.sh
# =====================================================================
set -euo pipefail

TEMPLATE=${TEMPLATE:-s2=0.15_s4=0.17}
NEW=${NEW:-s2=0.15_s4=0.15}
S2=${S2:-0.15}
S4=${S4:-0.15}
SUBMIT=${SUBMIT:-0}
JOBSCRIPT=${JOBSCRIPT:-run.sh}

echo "### template : $TEMPLATE"
echo "### new case : $NEW   (s2=$S2, s4=$S4)"
echo "### cwd      : $(pwd)"

[ -d "$TEMPLATE" ]            || { echo "!!! template dir $TEMPLATE not found"; exit 1; }
[ -f "$TEMPLATE/input.json" ] || { echo "!!! $TEMPLATE/input.json not found"; exit 1; }

if [ -d "$NEW" ]; then
  # second call after a dry run: submit the directory already built
  if [ "$SUBMIT" = "1" ] && [ ! -d "$NEW/snapshots" ]; then
    [ -f "$NEW/$JOBSCRIPT" ] || { echo "!!! $NEW/$JOBSCRIPT not found"; exit 1; }
    echo "### $NEW exists from the dry run and has no snapshots yet -- submitting it"
    ( cd "$NEW" && sbatch "$JOBSCRIPT" )
    exit 0
  fi
  echo "!!! $NEW already exists -- refusing to overwrite. Remove it first if you mean to rebuild."
  exit 1
fi

# --- copy only small top-level files (json, job script, helpers). Never
#     snapshots, checkpoints or logs from the template run.
mkdir "$NEW"
find "$TEMPLATE" -maxdepth 1 -type f -size -2M \
     ! -name '*.h5' ! -name '*.out' ! -name '*.err' ! -name '*.log' \
     -exec cp {} "$NEW"/ \;
echo "### copied files:"; ls -l "$NEW"

# --- edit s2, s4 in input.json; abort if those exact keys are missing
python3 - "$NEW/input.json" "$S2" "$S4" "$TEMPLATE" "$NEW" <<'EOF'
import json, sys
path, s2, s4, old, new = sys.argv[1], float(sys.argv[2]), float(sys.argv[3]), sys.argv[4], sys.argv[5]
d = json.load(open(path))
missing = [k for k in ("s2", "s4") if k not in d]
if missing:
    sys.exit(f"!!! keys {missing} not in input.json; keys present: {list(d)}. "
             "Edit by hand, nothing was changed.")
d["s2"], d["s4"] = s2, s4
# any string value that names the template directory (e.g. an output path)
renamed = []
for k, v in d.items():
    if isinstance(v, str) and old in v:
        d[k] = v.replace(old, new); renamed.append(k)
json.dump(d, open(path, "w"), indent=2)
print(f"### set s2={s2}, s4={s4}" + (f"; renamed paths in keys {renamed}" if renamed else ""))
# things that may need a per-case value -- print them for a human check
flag = [k for k in d if any(w in k.lower() for w in ("evp", "ic", "init", "sigma", "omega"))]
if flag:
    print("### CHECK BY HAND -- these keys may be specific to the template pair:")
    for k in flag: print(f"      {k} = {d[k]}")
EOF

echo
echo "### input.json diff (template -> new):"
diff <(python3 -m json.tool "$TEMPLATE/input.json") <(python3 -m json.tool "$NEW/input.json") || true

# --- other copied files that still mention the template name
echo
echo "### other files still naming the template (edit by hand if needed):"
grep -l -- "$TEMPLATE" "$NEW"/* 2>/dev/null | grep -v input.json || echo "      none"

# --- any EVP / eigenvector initial-condition files copied from the template
if ls "$NEW" | grep -qi 'evp\|eigen'; then
  echo
  echo "!!! An EVP / eigenvector file was copied from the template. The initial"
  echo "!!! condition depends on (s2,s4): recompute it for $NEW before submitting."
fi

if [ "$SUBMIT" = "1" ]; then
  [ -f "$NEW/$JOBSCRIPT" ] || { echo "!!! $NEW/$JOBSCRIPT not found"; exit 1; }
  echo; echo "### submitting $NEW/$JOBSCRIPT"
  ( cd "$NEW" && sbatch "$JOBSCRIPT" )
else
  echo; echo "### DRY RUN -- nothing submitted. Check the output above, then:"
  echo "      SUBMIT=1 bash $0"
fi
