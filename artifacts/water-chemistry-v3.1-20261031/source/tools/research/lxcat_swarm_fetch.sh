#!/usr/bin/env bash
# Fetch all e/<species> swarm/transport entries of one legacy LXCat database
# through the "Browse and download" form flow (nl.lxcat.net, data type "swrm")
# and store the raw Download.zip, the cache index page, the process-selection
# page and the request description in OUTDIR.
# Usage: lxcat_swarm_fetch.sh <db_id> <species_id> <outdir> <cookiejar>
# The cookie jar must already carry the accepted terms-of-use session.
set -euo pipefail
DB=$1; SPEC=$2; OUT=$3; CJ=$4
BASE=https://nl.lxcat.net
WORK=$(mktemp -d -p "${SCRATCH_DIR:-/tmp}")
trap 'rm -rf "$WORK"' EXIT
C="curl -sS -m 120 -L -c $CJ -b $CJ"
$C $BASE/data/ -o /dev/null
$C -X POST -d 'type=swrm&spec[]=elec' $BASE/data/set_type.php -o /dev/null
$C -X POST -d "db[]=$DB" $BASE/data/set_databases.php -o "$WORK/s3.html"
if grep -q 'name="spec\[\]"' "$WORK/s3.html"; then
  $C -X POST -d "spec[]=$SPEC" $BASE/data/set_specB.php -o "$WORK/s4.html"
else
  cp "$WORK/s3.html" "$WORK/s4.html"
fi
if grep -q 'name="gr\[\]"' "$WORK/s4.html"; then
  GR=$(grep -o '<input[^>]*name="gr\[\]"[^>]*>' "$WORK/s4.html" | sed -E 's/.*value="([^"]*)".*/-d gr[]=\1/' | tr '\n' ' ')
  $C -X POST $GR $BASE/data/set_groups.php -o "$WORK/s5.html"
else
  cp "$WORK/s4.html" "$WORK/s5.html"
fi
PR=$(grep -o '<input[^>]*name="proc\[\]"[^>]*>' "$WORK/s5.html" | sed -E 's/.*value="([^"]*)".*/-d proc[]=\1/' | tr '\n' ' ')
NPROC=$(grep -o '<input[^>]*name="proc\[\]"[^>]*>' "$WORK/s5.html" | wc -l)
URL=$($C -X POST $PR $BASE/data/set_processes.php -o "$WORK/s6.html" -w '%{url_effective}')
ID=${URL##*\?}
CACHE=""
for i in 1 2 3 4 5 6 7 8 9 10; do
  sleep 2
  CACHE=$($C "$BASE/data/output.php?$ID" -o "$WORK/s7.html" -w '%{url_effective}')
  case "$CACHE" in */cache/*) break;; esac
done
case "$CACHE" in */cache/*) ;; *) echo "no cache page for $DB" >&2; exit 2;; esac
T=$(date -u +%Y-%m-%dT%H:%M:%SZ)
$C "${CACHE}Download.zip" -o "$OUT/Download.zip"
cp "$WORK/s7.html" "$OUT/lxcat_output_page.html"
cp "$WORK/s5.html" "$OUT/lxcat_process_selection_page.html"
echo "$T" > "$OUT/retrieved_utc.txt"
{
  echo "LXCat legacy browse flow (nl.lxcat.net): type=swrm, spec[]=elec; db[]=$DB; spec[]=$SPEC; all data groups; all $NPROC processes selected."
  echo "Cache page: $CACHE"
  echo "Raw file: ${CACHE}Download.zip"
} > "$OUT/request.txt"
echo "$DB nproc=$NPROC cache=$CACHE"
