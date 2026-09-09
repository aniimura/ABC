"""Skeleton の主張が実際どんな形をしているかを数える。

問い: 主張は「A ならば B」の形か。
測り方: 各 theorem の statement 部分（`theorem` 行から `:=` の手前まで）を取り、
  ・仮説の本数（`(h... : ...)` 形の束縛）
  ・結論の先頭記号（∃ / ↔ / ∀ / ¬ / それ以外）
で分ける。★近似（構文的）であることは承知のうえで、比率を見る。
"""
import io
import os
import re
import collections

BS = chr(92)
ROOT = 'lean/ABC3/Skeleton'
DECL = re.compile(r'^(?:noncomputable\s+|private\s+)?(theorem|lemma)\s+'
                  r'([A-Za-z_][\w' + "'" + r'!?]*)')

files = {}
decls = []
for dp, _, fs in os.walk(ROOT):
    for f in fs:
        if not f.endswith('.lean'):
            continue
        p = os.path.join(dp, f).replace(BS, '/')
        t = io.open(p, encoding='utf-8').read()
        files[p] = t.split('\n')
        # 最初の namespace より後だけを宣言とみなす（#360）
        ns = 0
        for i, l in enumerate(files[p], 1):
            if l.startswith('namespace'):
                ns = i
                break
        for i, l in enumerate(files[p], 1):
            if i <= ns:
                continue
            m = DECL.match(l)
            if m and not m.group(2).endswith(('.src', '.needs')):
                decls.append((p, i, m.group(2)))


def statement_of(p, i):
    out = []
    for l in files[p][i - 1:i - 1 + 120]:
        j = l.find(':=')
        if j >= 0:
            out.append(l[:j])
            break
        out.append(l)
    return '\n'.join(out)


shape = collections.Counter()
hypcount = collections.Counter()
rows = []
for (p, i, n) in decls:
    st = statement_of(p, i)
    # 仮説とみなす束縛: `(名前 : ...)` で名前が h/H で始まるもの、または → を含む
    hyps = len(re.findall(r'\((h[A-Za-z0-9_]*|_h[A-Za-z0-9_]*)\s*:', st))
    arrows = st.count('→')
    # 結論部（最後の `:` の後ろ）をざっくり取る
    tail = st.rsplit(':', 1)[-1]
    if '↔' in st:
        s = 'iff（同値）'
    elif re.search(r'∃', tail):
        s = 'exists（存在）'
    elif re.search(r'¬', tail):
        s = 'not（否定）'
    elif arrows > 0 or hyps > 0:
        s = 'implication（A ならば B）'
    else:
        s = 'bare（仮説なしの等式・不等式など）'
    shape[s] += 1
    hypcount[min(hyps, 5)] += 1
    rows.append((n, s, hyps, p.replace('lean/ABC3/Skeleton/', '')))

print('Skeleton の theorem/lemma: %d 件' % len(decls))
print()
print('■ 形の分布')
for k, v in shape.most_common():
    print('   %-28s %4d  (%.0f%%)' % (k, v, 100.0 * v / len(decls)))
print()
print('■ 仮説（h… 束縛）の本数')
for k in sorted(hypcount):
    lab = '%d 本' % k if k < 5 else '5 本以上'
    print('   %-10s %4d' % (lab, hypcount[k]))
print()
print('■ bare（仮説なし）の例 10 件')
c = 0
for n, s, h, f in rows:
    if s == 'bare（仮説なしの等式・不等式など）' and c < 10:
        print('   %-44s %s' % (n, f))
        c += 1
