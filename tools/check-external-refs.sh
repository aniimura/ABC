#!/usr/bin/env bash
# ★external/_refs が lake-manifest.json と一致しているか確認する(数秒)。
# ズレていたら tools/sync-external-refs.sh を実行すること。
set -euo pipefail
cd "$(dirname "$0")/.."
'C:/Users/Aruta/miniforge3/envs/py311env/python.exe' - <<'PY'
import json, io, os, sys
man = json.load(io.open('lean/lake-manifest.json', encoding='utf-8'))
cur = {p['name']: p.get('rev', '') for p in man.get('packages', [])}
sp = 'external/_refs/STAMP.json'
if not os.path.exists(sp):
    print('NG external/_refs が無い —— tools/sync-external-refs.sh を実行すること')
    sys.exit(1)
old = json.load(io.open(sp, encoding='utf-8'))['revs']
bad = [k for k, v in old.items() if cur.get(k) != v]
if bad:
    print('NG rev がズレている: ' + ', '.join(bad))
    print('   tools/sync-external-refs.sh を実行すること')
    sys.exit(1)
print('ok external/_refs は lake-manifest.json と一致(' +
      ', '.join('%s=%s' % (k, v[:12]) for k, v in sorted(old.items())) + ')')
PY
