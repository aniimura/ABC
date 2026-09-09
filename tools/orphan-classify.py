"""孤児（誰にも参照されていない Skeleton の宣言）を分類する。

goal-chain.json の skeletonOrphans を、木の実ファイルと突き合わせて
  ・種別（theorem / def / structure …）
  ・`.src` を持つか（= 原典の項目そのものの主張として登記されているか）
  ・どのファイルに居るか
で分ける。★削除してよいかはこの分類で決まる。
"""
import json
import io
import os
import re
import collections

BS = chr(92)

d = json.load(io.open('lean/.cache/goal-chain.json', encoding='utf-8'))
orph = d['skeletonOrphans']

DECL = re.compile(r'^(?:noncomputable\s+)?(theorem|lemma|def|abbrev|structure)\s+'
                  r'([A-Za-z_][\w.' + "'" + r'!?]*)')

src = {}
for dp, _, fs in os.walk('lean/ABC3/Skeleton'):
    for f in fs:
        if not f.endswith('.lean'):
            continue
        p = os.path.join(dp, f).replace(BS, '/')
        t = io.open(p, encoding='utf-8').read()
        for i, l in enumerate(t.split('\n'), 1):
            m = DECL.match(l)
            if m:
                src.setdefault(m.group(2), []).append((p, i, m.group(1)))


def short(n):
    return n.split('.')[-1]


cnt = collections.Counter()
rows = []
for n in orph:
    s = short(n)
    loc = src.get(s, [])
    kind = loc[0][2] if loc else '?'
    f = loc[0][0].replace('lean/ABC3/Skeleton/', '') if loc else '?'
    hassrc = (s + '.src') in src
    rows.append((n, kind, hassrc, f))
    cnt[(kind, hassrc)] += 1

print('孤児 %d 件 —— 種別 x .src の有無' % len(orph))
for k, v in sorted(cnt.items(), key=lambda x: -x[1]):
    print('   %-10s .src=%-5s  %d' % (k[0], k[1], v))
print()
print('★.src を持つ theorem/lemma（原典の項目そのものの主張。削除は形式化の破棄になる）')
for n, k, h, f in rows:
    if h and k in ('theorem', 'lemma'):
        print('   %-55s %s' % (n, f))
print()
print('★.src を持たない theorem/lemma（補助補題。消費者が無ければ削除候補）')
for n, k, h, f in rows:
    if (not h) and k in ('theorem', 'lemma'):
        print('   %-55s %s' % (n, f))
print()
print('★def/abbrev/structure（述語・データ。削除候補だが他が型で使う可能性あり）')
for n, k, h, f in rows:
    if k in ('def', 'abbrev', 'structure'):
        print('   %-55s %s  .src=%s' % (n, f, h))
print()
print('★木に見つからなかったもの（名前空間や where の中）')
for n, k, h, f in rows:
    if k == '?':
        print('   %s' % n)
