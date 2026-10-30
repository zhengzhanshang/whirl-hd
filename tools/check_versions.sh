#!/usr/bin/env bash
# =====================================================================
# check_versions.sh -- for the closed-form recovery method paper:
#
#   1. REPORT: find every copy of the paper's scripts on THIS machine,
#      group identical copies by checksum, and test each distinct version
#      for the features the paper's version must have (OK / OLD).
#
#   2. COLLECT: copy the archive-ready files into
#         <parent of cwd>/scripts_for_papers_data_avail/
#      - analysis / plotting / post-processing scripts: the newest copy
#        that passes every feature check (newest by modification time if
#        a script has no feature checks);
#      - Dedalus scripts, run.sh and input.json: taken from the paper's own
#        case folders, separated by G (the copy in a case folder is the one
#        that produced that run's data, so it is the one to archive);
#      - one input_template.json + pairs.txt per G when all input.json of
#        that G differ only in s2 and s4, otherwise one folder per pair;
#      - MANIFEST.txt with the checksum and source path of every copy.
#
# Run once on each machine (IGAM12WS holds the G=0.25 runs, GSC the
# G=4.5 runs). Each machine builds its own folder; merge the two by hand
# and use MANIFEST.txt to resolve any script found on both.
#
#   bash check_versions.sh                  # report + collect
#   COLLECT=0 bash check_versions.sh        # report only
#   ROOTS="/home/zshang/dedalus3" bash check_versions.sh
#   CASE_DIRS="/path/to/sweep1 /path/to/sweep2" bash check_versions.sh
#   PAIR_DIRS=1 bash check_versions.sh      # always one folder per pair
#
# PINNED versions: a script whose md5 is listed in PINS is the exact copy
# that produced the paper's numbers. The report marks it PINNED and the
# collect step takes it in preference to a newer copy.
#   lsq_from_data.py 176d096ac17a... = Table 1, Sect. 4.1 and the Sect. 4
#   verification reruns on IGAM12WS, 2026-10-19.
#
# The destination folder is never overwritten: if it exists, the collect
# step stops. Remove or rename it first to rebuild.
# =====================================================================
ROOTS=${ROOTS:-"$HOME /cl_tmp/$USER"}
OUT=${OUT:-versions_$(hostname -s).txt}
COLLECT=${COLLECT:-1}
PAIR_DIRS=${PAIR_DIRS:-0}
DEST=${DEST:-"$(dirname "$PWD")/scripts_for_papers_data_avail"}
# the paper's sweep folders; only those that exist on this machine are used
PINS=${PINS:-"lsq_from_data.py=176d096ac17a3371dc06f164e662de14"}
export PINS
CASE_DIRS=${CASE_DIRS:-"/home/zshang/dedalus3/HD_SW_alpha_spherical_paper/only_h_small_G=0.25_paper/vary_s2s4_G=0.25_nu4-h=1e-4_paper/sweep_s2s4_pairs /cl_tmp/shang/paper-redo_G=4.5/m=3-15/m=5_s2s4_sweep"}

SCRIPTS=(spherical-2d-hd-s2s4.py sw_model.py evp_spherical_hd.py
         postrun_analyse.py postrun_analyze.py measure_mode.py fit_growth.py
         lsq_from_data.py los_inversion.py lsq_hfree.py make_figures.py
         select_settings.py)
DRIVERS=(run_table1_G025.sh run_pole_mask_G025.sh run_verify_numbers_G025.sh
         run_sec41_mask.sh run_sec44_cadence.sh
         run_sec5_los.sh run_gdep_G4.5_m5.sh run_hfree_los.sh
         make_case_s015_s015.sh check_versions.sh)

# "file|marker|meaning": a marker the paper's version MUST contain.
# A leading ! means the marker must be ABSENT.
MARKERS=(
 "lsq_from_data.py|--ckpt|reads inverted fields (Sect. 5)"
 "lsq_from_data.py|eq-subset|equation subsets (Sect. 4.5)"
 "lsq_from_data.py|--dropout|random dropout (Sect. 4.3)"
 "lsq_from_data.py|--gap|contiguous gap (Sect. 4.3)"
 "lsq_from_data.py|tstride|cadence test (Sect. 4.4)"
 "lsq_from_data.py|limb|centre-to-limb weighting (Sect. 4.2)"
 "lsq_from_data.py|mask-poles|polar masking (Sect. 4.1)"
 "lsq_from_data.py|nu4_h|hyperdiffusion in eq3 (Sect. 2.2)"
 "lsq_from_data.py|Fornberg|non-uniform time weights (Sect. 2.4)"
 "lsq_from_data.py|TRUTH-FREE|truth-free selector (Sect. 2.5)"
 "los_inversion.py|amp|amplitude-threshold lambda fit (Sect. 5.3)"
 "los_inversion.py|pinn_fields.h5|writes checkpoints lsq can read"
 "lsq_hfree.py|--band-deg|band-limited derivatives (Sect. 2.4)"
 "lsq_hfree.py|band-trim|band-edge trimming (Sect. 2.4)"
 "lsq_hfree.py|--no-visc|hyperdiffusion switch (Sect. 2.4)"
 "lsq_hfree.py|nu4=0 limit|selftest covers nu4 terms (App. B)"
 "make_figures.py|eq1+2+3|Fig. 2 legend fixed"
 "make_figures.py|!eq1+2+5|old Fig. 2 legend gone"
 "make_figures.py|!def fig4|PINN figure removed"
 "make_figures.py|PAIRS|three-pair Fig. 1"
 "make_figures.py|hspace=0.32|Fig. 3 panel spacing fixed"
 "evp_spherical_hd.py|--write-ic|writes the initial condition"
 "evp_spherical_hd.py|--pick|selects the Rossby mode"
 "spherical-2d-hd-s2s4.py|evp_ic|starts from the EVP eigenvector"
 "spherical-2d-hd-s2s4.py|nu4_h|hyperdiffusion on h"
 "postrun_analyse.py|sigma_spread|fixed clean-mode verdict"
)

TMP=$(mktemp -d)
trap 'rm -rf "$TMP"' EXIT

# find copies, never descending into a previously built archive folder
findcopies() { find $ROOTS -path "$DEST" -prune -o -name "$1" -type f -print 2>/dev/null | sort; }

exec > >(tee "$OUT") 2>&1
echo "# check_versions  host=$(hostname)  user=$USER  $(date -Is)"
echo "# roots: $ROOTS"

# ------------------------------------------------------------- 1. REPORT
for S in "${SCRIPTS[@]}" "${DRIVERS[@]}"; do
  mapfile -t FILES < <(findcopies "$S")
  echo
  echo "====================================================================="
  echo "== $S   (${#FILES[@]} copies)"
  echo "====================================================================="
  [ ${#FILES[@]} -eq 0 ] && { echo "   none found"; continue; }

  echo "-- every copy, newest first: modified | size | md5 | path"
  for f in "${FILES[@]}"; do
    printf "%s|%s|%s|%s\n" "$(date -r "$f" '+%Y-%m-%d %H:%M')" \
      "$(stat -c %s "$f")" "$(md5sum "$f" | cut -c1-12)" "$f"
  done | sort -r | awk -F'|' '{printf "   %-16s %9s  %s  %s\n", $1, $2, $3, $4}'

  echo "-- distinct versions (copies per md5)"
  for f in "${FILES[@]}"; do md5sum "$f"; done | cut -c1-12 | sort | uniq -c \
    | sort -rn | sed 's/^/   /'

  # one representative path per distinct version
  mapfile -t REPS < <(for f in "${FILES[@]}"; do echo "$(md5sum "$f" | cut -c1-12) $f"; done \
                      | sort -u -k1,1 | cut -d' ' -f2-)
  HAVE=0
  for M in "${MARKERS[@]}"; do [[ "$M" == "$S|"* ]] && HAVE=1; done
  if [ $HAVE -eq 0 ]; then
    # no feature checks: every version is a candidate (newest wins later)
    printf "%s\n" "${REPS[@]}" > "$TMP/cand_$S"
    echo "-- no feature checks for this file; the newest copy will be collected"
    continue
  fi
  echo "-- feature check per version (OK = matches the paper's version)"
  : > "$TMP/cand_$S"
  for f in "${REPS[@]}"; do
    miss=0; line=""
    for M in "${MARKERS[@]}"; do
      [[ "$M" == "$S|"* ]] || continue
      IFS='|' read -r _ pat what <<< "$M"
      if [[ "$pat" == !* ]]; then
        grep -qF -- "${pat:1}" "$f" && { miss=1; line+=" [still contains '${pat:1}', should not: $what]"; }
      else
        grep -qF -- "$pat" "$f" || { miss=1; line+=" [missing '$pat': $what]"; }
      fi
    done
    pin=""
    [[ " $PINS " == *" $S=$(md5sum "$f" | cut -d' ' -f1)"* ]] && pin="   <-- PINNED: produced the paper's numbers"
    if [ $miss -eq 0 ]; then
      echo "   OK   $(md5sum "$f" | cut -c1-12)  $f$pin"
      echo "$f" >> "$TMP/cand_$S"
    else
      echo "   OLD  $(md5sum "$f" | cut -c1-12)  $f$pin"; echo "        $line"
      # the pinned copy produced the paper's numbers: collect it even so
      [ -n "$pin" ] && { echo "        (pinned: collected anyway -- the markers list may be stricter than needed)"; echo "$f" >> "$TMP/cand_$S"; }
    fi
  done
done

echo
echo "# report done. Compare this file with the one from the other machine:"
echo "#   identical md5 = identical file; archive one copy of each OK version."

# ------------------------------------------------------------- 2. COLLECT
[ "$COLLECT" = "1" ] || { echo "# COLLECT=0 -- nothing copied."; exit 0; }

echo
echo "#####################################################################"
echo "# COLLECT into $DEST"
echo "#####################################################################"
if [ -e "$DEST" ]; then
  echo "!!! $DEST already exists -- not touching it. Remove or rename it to rebuild."
  exit 1
fi

python3 - "$DEST" "$TMP" "$PAIR_DIRS" "$(hostname)" $CASE_DIRS <<'PYEOF'
import glob, hashlib, json, math, os, shutil, sys, time

dest, tmp, pair_dirs, host = sys.argv[1], sys.argv[2], sys.argv[3] == "1", sys.argv[4]
case_roots = [d for d in sys.argv[5:] if os.path.isdir(d)]

GROUPS = {
    "analysis":    ["lsq_from_data.py", "los_inversion.py", "lsq_hfree.py"],
    "figures":     ["make_figures.py"],
    "postprocess": ["postrun_analyse.py", "measure_mode.py", "fit_growth.py"],
    "drivers":     ["run_table1_G025.sh", "run_pole_mask_G025.sh",
                    "run_verify_numbers_G025.sh", "run_sec41_mask.sh", "run_sec44_cadence.sh",
                    "run_sec5_los.sh", "run_gdep_G4.5_m5.sh", "run_hfree_los.sh",
                    "make_case_s015_s015.sh", "check_versions.sh"],
}
SIM = ["spherical-2d-hd-s2s4.py", "sw_model.py", "evp_spherical_hd.py", "run.sh"]
PAPER_PAIRS = {
    0.25: [(0.29, 0.03), (0.25, 0.19), (0.19, 0.23), (0.15, 0.17), (0.15, 0.15), (0.03, 0.33)],
    4.5:  [(0.29, 0.03), (0.25, 0.19), (0.19, 0.23), (0.15, 0.17), (0.15, 0.15), (0.16, 0.12),
           (0.03, 0.33)],
}
M_PAPER = 5
PINS = dict(x.split("=", 1) for x in os.environ.get("PINS", "").split() if "=" in x)

manifest, warnings = [], []
def md5(p):
    return hashlib.md5(open(p, "rb").read()).hexdigest()
def mtime(p):
    return time.strftime("%Y-%m-%d %H:%M", time.localtime(os.path.getmtime(p)))
def copy(src, rel, note=""):
    out = os.path.join(dest, rel)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    shutil.copy2(src, out)
    manifest.append((rel, md5(out), src, mtime(src), note))

os.makedirs(dest)

# ---- scripts that are not tied to a run: newest OK version
for sub, names in GROUPS.items():
    for name in names:
        cf = os.path.join(tmp, "cand_" + name)
        cands = [l.strip() for l in open(cf)] if os.path.exists(cf) else []
        cands = [c for c in cands if c and os.path.isfile(c)]
        if not cands:
            warnings.append(f"{name}: no OK copy found on {host} -- NOT collected")
            continue
        cands.sort(key=os.path.getmtime, reverse=True)
        note = ""
        pinned = PINS.get(name)
        if pinned:
            hit = [c for c in cands if md5(c) == pinned]
            if hit:
                cands = hit + [c for c in cands if c not in hit]
                note = "PINNED version (produced the paper's numbers). "
            else:
                warnings.append(f"{name}: pinned md5 {pinned[:12]} NOT found on {host}; newest OK copy taken")
        if len(cands) > 1 and not (pinned and note):
            note += (f"{len(cands)} distinct OK versions; newest taken. Others: "
                    + "; ".join(cands[1:]) + "  -- check with diff")
            warnings.append(f"{name}: {note}")
        copy(cands[0], os.path.join(sub, name), note)

# ---- the paper's simulation cases, by G
def gkey(G):
    for g in PAPER_PAIRS:
        if math.isclose(G, g, rel_tol=1e-9, abs_tol=1e-12):
            return g
    return None

cases = {}   # (G, s2, s4) -> case dir
for root in case_roots:
    for inp in sorted(glob.glob(os.path.join(root, "*", "input.json"))):
        try:
            d = json.load(open(inp))
            G, m, s2, s4 = float(d["G"]), int(d["m"]), float(d["s2"]), float(d["s4"])
        except Exception as ex:
            warnings.append(f"{inp}: unreadable or missing G/m/s2/s4 ({ex}) -- skipped")
            continue
        g = gkey(G)
        if g is None or m != M_PAPER:
            continue
        pair = next((p for p in PAPER_PAIRS[g]
                     if math.isclose(p[0], s2) and math.isclose(p[1], s4)), None)
        if pair is None:
            continue
        key = (g, pair)
        cdir = os.path.dirname(inp)
        if key in cases:
            if md5(inp) != md5(os.path.join(cases[key], "input.json")):
                warnings.append(f"G={g} {pair}: two case folders with DIFFERENT input.json: "
                                f"{cases[key]} (kept) and {cdir} (ignored)")
            continue
        cases[key] = cdir

for g, pairs in PAPER_PAIRS.items():
    found = [(p, cases[(g, p)]) for p in pairs if (g, p) in cases]
    if not found:
        continue
    gdir = f"G{g:g}_m{M_PAPER}"
    missing = [p for p in pairs if (g, p) not in cases]
    if missing:
        warnings.append(f"{gdir}: pairs not found on {host}: {missing}")

    # simulation scripts: one shared copy if identical across the pairs
    for f in SIM:
        have = [(p, os.path.join(c, f)) for p, c in found if os.path.isfile(os.path.join(c, f))]
        if not have:
            warnings.append(f"{gdir}: {f} not found in any case folder")
            continue
        sums = {md5(path) for _, path in have}
        if len(sums) == 1:
            copy(have[0][1], os.path.join(gdir, "simulation", f),
                 f"identical in all {len(have)} case folders")
        else:
            warnings.append(f"{gdir}: {f} differs between pairs -- one copy per pair")
            for p, path in have:
                copy(path, os.path.join(gdir, f"s2={p[0]:g}_s4={p[1]:g}", f),
                     "differs between pairs")

    # input.json: template if they differ only in s2 and s4
    def norm(path, cdir):
        d = json.load(open(path))
        name = os.path.basename(cdir)
        out = {}
        for k, v in d.items():
            if k in ("s2", "s4"):
                continue
            out[k] = v.replace(name, "<PAIR>") if isinstance(v, str) else v
        return json.dumps(out, sort_keys=True)
    normed = {norm(os.path.join(c, "input.json"), c) for _, c in found}
    if len(normed) == 1 and not pair_dirs:
        p0, c0 = found[0]
        tmpl = json.load(open(os.path.join(c0, "input.json")))
        tmpl["s2"] = "SET_PER_PAIR"; tmpl["s4"] = "SET_PER_PAIR"
        name = os.path.basename(c0)
        tmpl = {k: (v.replace(name, "<PAIR>") if isinstance(v, str) else v)
                for k, v in tmpl.items()}
        out = os.path.join(dest, gdir, "input_template.json")
        os.makedirs(os.path.dirname(out), exist_ok=True)
        json.dump(tmpl, open(out, "w"), indent=2)
        manifest.append((os.path.join(gdir, "input_template.json"), md5(out),
                         os.path.join(c0, "input.json"), mtime(os.path.join(c0, "input.json")),
                         "template; s2, s4 set per pair from pairs.txt; <PAIR> = case folder name"))
        with open(os.path.join(dest, gdir, "pairs.txt"), "w") as fh:
            fh.write(f"# {gdir}: every input.json is input_template.json with these s2, s4\n")
            fh.write("# s2      s4      case folder\n")
            for p, c in found:
                fh.write(f"{p[0]:<8g}{p[1]:<8g}{os.path.basename(c)}\n")
    else:
        if len(normed) > 1:
            dicts = [json.loads(norm(os.path.join(c, "input.json"), c)) for _, c in found]
            keys = sorted({k for d in dicts for k in d}
                          - {k for k in dicts[0] if all(d.get(k) == dicts[0][k] for d in dicts)})
            warnings.append(f"{gdir}: input.json differ beyond s2/s4 (keys: {', '.join(keys)})"
                            " -- one folder per pair")
        for p, c in found:
            copy(os.path.join(c, "input.json"),
                 os.path.join(gdir, f"s2={p[0]:g}_s4={p[1]:g}", "input.json"))

# ---- manifest
with open(os.path.join(dest, "MANIFEST.txt"), "w") as fh:
    fh.write(f"# scripts_for_papers_data_avail -- built on {host} {time.strftime('%Y-%m-%d %H:%M')}\n")
    fh.write("# archived path | md5 | copied from | source modified | note\n\n")
    for rel, h, src, mt, note in sorted(manifest):
        fh.write(f"{rel}\n    md5 {h}\n    from {src}  ({mt})\n")
        if note:
            fh.write(f"    note: {note}\n")
    if warnings:
        fh.write("\n# WARNINGS -- check these before uploading\n")
        for w in warnings:
            fh.write(f"!!! {w}\n")

print(f"copied {len(manifest)} files into {dest}")
for w in warnings:
    print("!!!", w)
print(f"see {os.path.join(dest, 'MANIFEST.txt')}")
PYEOF
echo "# collect done."
