#!/usr/bin/env bash
# ★既存数学を「既定の grep 範囲」に置くための同期スクリプト(2026-08-19)。
#
# なぜ要るか:
#   mathlib は lean/.lake/packages/ にあるが、lean/.lake/ が .gitignore で
#   ディレクトリごと外れているため、**既定の grep に入らない**。
#   gitignore の意味論で「親が除外されていると中身は再包含できない」ので
#   .ignore の否定では戻せない(2026-08-19 実測)。
#   ジャンクション(mklink /J)も駄目 —— ripgrep は --follow なしで辿らない(実測)。
#
#   そこで **.lean だけを external/ へ複製**する。external/ は .ignore で
#   打ち消してあるので既定検索に入る。git は .gitignore で無視したままなので
#   再配布もしないし、ABC3b の `git add -A` にも巻き込まれない。
#
# 陳腐化への備え:
#   lake-manifest.json の rev を STAMP に書き込む。ズレたら再実行すること。
#   確認は tools/check-external-refs.sh(rev 比較のみ、数秒)。
set -euo pipefail
cd "$(dirname "$0")/.."
SRC=lean/.lake/packages
DST=external/_refs
mkdir -p "$DST"
for p in mathlib batteries PrimeNumberTheoremAnd; do
  case "$p" in
    mathlib) sub=Mathlib ;;
    batteries) sub=Batteries ;;
    PrimeNumberTheoremAnd) sub=PrimeNumberTheoremAnd ;;
  esac
  if [ -d "$SRC/$p/$sub" ]; then
    rm -rf "$DST/$p"
    mkdir -p "$DST/$p"
    (cd "$SRC/$p" && find "$sub" -name "*.lean" -print0) \
      | (cd "$SRC/$p" && tar --null -cf - --files-from=-) \
      | (cd "$DST/$p" && tar -xf -)
    echo "  $p/$sub -> $DST/$p ($(find "$DST/$p" -name '*.lean' | wc -l) files)"
  fi
done
python_exe='C:/Users/Aruta/miniforge3/envs/py311env/python.exe'
"$python_exe" - <<'PY'
import json, io, datetime
d = json.load(io.open('lean/lake-manifest.json', encoding='utf-8'))
rows = [(p['name'], p.get('rev', '')) for p in d.get('packages', [])
        if p['name'] in ('mathlib', 'batteries', 'PrimeNumberTheoremAnd')]
out = io.open('external/_refs/STAMP.json', 'w', encoding='utf-8', newline='')
out.write(json.dumps({'syncedFrom': 'lean/lake-manifest.json',
                      'revs': dict(rows)}, ensure_ascii=False, indent=1) + '\n')
out.close()
print('  STAMP.json:', dict((k, v[:12]) for k, v in rows))
PY
