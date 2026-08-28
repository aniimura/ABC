import io, re, json, os

ROOT = r'D:\Math_ABC3\lean\ABC3'

# item -> (skeleton file, theorem name)
TARGETS = [
    ("Proposition 1.4",  "Skeleton/GenEll/Section1.lean", "prop_1_4"),
    ("Remark 1.4.1",     "Skeleton/GenEll/Section1.lean", "remark_1_4_1"),
    ("Proposition 1.7",  "Skeleton/GenEll/Section1.lean", "prop_1_7"),
    ("Theorem 2.1",      "Skeleton/GenEll/Section2.lean", "theorem_2_1"),
    ("Lemma 3.2",        "Skeleton/GenEll/Section3.lean", "lemma_3_2"),
    ("Proposition 3.4",  "Skeleton/GenEll/Section3.lean", "prop_3_4"),
    ("Lemma 3.5",        "Skeleton/GenEll/Section3.lean", "lemma_3_5"),
    ("Lemma 3.7",        "Skeleton/GenEll/Section3.lean", "lemma_3_7"),
    ("Theorem 3.8",      "Skeleton/GenEll/GaloisImage.lean", "theorem_3_8"),
    ("Corollary 4.3",    "Skeleton/GenEll/Section4.lean", "cor_4_3"),
    ("Corollary 4.4",    "Skeleton/GenEll/Section4.lean", "cor_4_4"),
]

def body_of(path, name):
    src = io.open(os.path.join(ROOT, path), encoding='utf-8').read()
    m = re.search(r'^theorem\s+' + re.escape(name) + r'\b', src, re.M)
    if m is None:
        return None
    start = m.start()
    # next top-level declaration or section marker
    nxt = re.search(r'^(theorem|def|noncomputable def|/-!|end )', src[m.end():], re.M)
    end = m.end() + (nxt.start() if nxt else len(src) - m.end())
    return src[start:end]

rows = []
for item, path, name in TARGETS:
    b = body_of(path, name)
    if b is None:
        rows.append((item, path, name, None, None, None, None))
        continue
    has_sorry = 'sorry' in b
    # interface projections consumed: D.foo / S.foo occurrences
    projs = sorted(set(re.findall(r'\bD\.([a-zA-Z_][a-zA-Z_0-9\']*)', b)))
    # calls into Found/
    found = sorted(set(re.findall(r'ABC3\.Found\.GenEll\.([a-zA-Z_][a-zA-Z_0-9\'_]*)', b)))
    # calls into other Skeleton theorems
    skel = sorted(set(re.findall(r'\b(theorem_3_8|prop_3_4|lemma_3_[0-9]|cor_4_[0-9]|lemma_4_[0-9])\b', b)))
    skel = [t for t in skel if t != name]
    rows.append((item, path, name, has_sorry, projs, found, skel))

out = []
out.append("# GenEll 残り 11 項目の実測（2026-08-29）\n")
out.append("| 項目 | Skeleton | sorry | Interface 射影 | Found 呼び出し | 他の Skeleton 定理 |")
out.append("|---|---|---|---|---|---|")
for item, path, name, has_sorry, projs, found, skel in rows:
    if has_sorry is None:
        out.append(f"| {item} | {name}（未検出） | — | — | — | — |")
        continue
    out.append("| {} | `{}` | {} | {} | {} | {} |".format(
        item, name, "あり" if has_sorry else "なし",
        len(projs) if projs else 0,
        len(found) if found else 0,
        ", ".join(f"`{t}`" for t in skel) if skel else "—"))
out.append("")
out.append("## 内訳\n")
for item, path, name, has_sorry, projs, found, skel in rows:
    if has_sorry is None:
        continue
    out.append(f"### {item}（`{name}`、sorry {'あり' if has_sorry else 'なし'}）\n")
    out.append("* Interface 射影: " + (", ".join(f"`D.{p}`" for p in projs) if projs else "なし"))
    out.append("* Found 呼び出し: " + (", ".join(f"`{p}`" for p in found) if found else "なし"))
    out.append("* 他の Skeleton 定理: " + (", ".join(f"`{t}`" for t in skel) if skel else "なし"))
    out.append("")

io.open(r'C:\Users\Aruta\.claude\jobs\663400ce\tmp\remaining-measure.md', 'w',
         encoding='utf-8', newline='\n').write('\n'.join(out))
print('\n'.join(out[:20]))
