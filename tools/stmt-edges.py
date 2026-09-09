"""Skeleton の主張どうしが「証明を書く前に」型で繋がっているかを数える。

問い: 依存関係は証明まで形式化しないと判明しないのか。
測り方: 各 theorem の **statement 部分だけ**（`theorem` 行から `:=` の直前まで）を取り、
そこに別の Skeleton 主張の名前が現れるかを数える。
  ・現れる  → 証明を書かなくても依存が型で見える
  ・現れない → その依存は証明の中にしか無い（＝証明を書くまで機械には見えない）
"""
import io
import os
import re
import collections

BS = chr(92)
ROOT = 'lean/ABC3/Skeleton'
DECL = re.compile(r'^(?:noncomputable\s+|private\s+)?(theorem|lemma|def|abbrev|structure)\s+'
                  r'([A-Za-z_][\w' + "'" + r'!?]*)')

files = {}
decls = []          # (file, line, kind, name)
for dp, _, fs in os.walk(ROOT):
    for f in fs:
        if not f.endswith('.lean'):
            continue
        p = os.path.join(dp, f).replace(BS, '/')
        t = io.open(p, encoding='utf-8').read()
        files[p] = t.split('\n')
        for i, l in enumerate(files[p], 1):
            m = DECL.match(l)
            if m:
                decls.append((p, i, m.group(1), m.group(2)))

has_src = set()
for p, i, k, n in decls:
    if n.endswith('.src'):
        has_src.add(n[:-4])
# `.src` は `def foo.src` の形なので上の DECL では名前に `.` が入らない。別に拾う。
SRC = re.compile(r'^\s*def\s+([A-Za-z_][\w.' + "'" + r'!?]*)\.src\b')
for p, lines in files.items():
    for l in lines:
        m = SRC.match(l)
        if m:
            has_src.add(m.group(1))

# 主張とみなすもの: theorem/lemma で `.src` を持つ（＝原典の項目）
claims = [(p, i, n) for (p, i, k, n) in decls
          if k in ('theorem', 'lemma') and n in has_src and not n.endswith('.src')]
claim_names = set(n for (_, _, n) in claims)


def statement_of(p, i):
    """`theorem` 行から `:=` / `:= by` の直前までを返す（statement 部分だけ）。"""
    lines = files[p]
    out = []
    for l in lines[i - 1:i - 1 + 200]:
        j = l.find(':=')
        if j >= 0:
            out.append(l[:j])
            break
        out.append(l)
        if re.match(r'^(theorem|lemma|def|abbrev|structure|end|namespace)\b', l) and len(out) > 1:
            break
    return '\n'.join(out)


WORD = re.compile(r"[A-Za-z_][\w'!?]*")
linked = 0
edges = collections.Counter()
examples = []
for (p, i, n) in claims:
    st = statement_of(p, i)
    hits = set(w for w in WORD.findall(st) if w in claim_names and w != n)
    if hits:
        linked += 1
        edges[p.replace('lean/ABC3/Skeleton/', '')] += len(hits)
        if len(examples) < 12:
            examples.append((n, sorted(hits)[:3], p.replace('lean/ABC3/Skeleton/', '')))

print('Skeleton の主張（theorem/lemma で `.src` あり）: %d 件' % len(claims))
print('★そのうち statement の中に別の主張の名前が現れるもの: %d 件' % linked)
print('★現れないもの（依存が証明の中にしか無い）: %d 件' % (len(claims) - linked))
print()
if examples:
    print('statement で繋がっている例:')
    for n, h, f in examples:
        print('   %-42s -> %s   [%s]' % (n, ', '.join(h), f))
print()
print('ファイル別の statement 辺の数（上位 10）:')
for f, c in edges.most_common(10):
    print('   %-52s %d' % (f, c))

# ★構造上の理由を測る: 主張どうしが「共通の Data を取る兄弟」になっていないか
IFACE = re.compile(r'\b[A-Z][\w]*(?:Data|Setup|Filtration|Cardinality|Correspondence)\b')
sib = sum(1 for (p, i, n) in claims if IFACE.search(statement_of(p, i)))
print()
print('★statement が Interface のデータ構造（…Data / …Setup 等）を取る主張: %d / %d'
      % (sib, len(claims)))
