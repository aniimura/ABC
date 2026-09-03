import fs from 'fs';
const p = 'D:/Math_ABC3/ResearchPaper/blocked-leaves.json';
const j = JSON.parse(fs.readFileSync(p, 'utf8'));
const key = '[GenEll] Lemma 3.5(Vélu の商の半安定性——判別式の恒等式と p ∣ l の形式群)';
const e = j.blocked.find(b => b.key === key);
if (!e) { console.log('not found'); process.exit(1); }
e.identityClosed2026_09_02 =
  '★★★★★★★★★★★★★★★★★★★★**判別式の恒等式は閉じた（2026-09-02、第 1397-1402）**。'
  + '☆`exactMissingLemma` の（1）が無くなった——`Skeleton/GenEll/VeluDiscIdentity.lean` は '
  + '`sorry` 0 になり、`disc_pow_eq_veluQuot_mul`（数体版）と '
  + '`semistableAt_veluQuot_good`（良い素点、`p ∤ l`）が**無条件**になった。'
  + '★測った 3 つの道（σ 函数／終結式／剰余体の Vélu）はどれも要らず、'
  + '**Liouville（第 598、在庫）と `R` の偶関数性だけ**で出た'
  + '（`mathlib-gap.json` の `veluDiscIdentityClosed20260902` に全文）。'
  + '☆残るのは（a）組み立て（悪い素点の局所パッケージ、良い素点の極小モデルへの正規化）と'
  + '（b）**形式群の鋭い捩れ限界** `(l−1)·k ≤ v_p(l)`（`p ∣ l` の場合）の 2 つである。';
fs.writeFileSync(p, JSON.stringify(j, null, 2) + '\n', 'utf8');
console.log('recorded blocked-leaves');
