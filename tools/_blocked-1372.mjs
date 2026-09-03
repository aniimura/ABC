import fs from 'fs';
const p = 'D:/Math_ABC3/ResearchPaper/blocked-leaves.json';
const j = JSON.parse(fs.readFileSync(p, 'utf8'));
const e = j.blocked.find(b => b.key.startsWith('[GenEll] Theorem 3.8(悪い素点の惰性'));
if (!e) { console.log('not found'); process.exit(1); }
e.resolved2026_09_02 = '★★★★★★★★**この節点の mathlib 欠落は塞がった（2026-09-02、第 1365-1372）**。'
  + '☆第 1356 で特定した primitive「完備なイデアルに沿う冪等元の持ち上げ」を第 1357 で作り、'
  + '第 1358-1364 で 8 段を積み、第 1365（整拡大では `m_R·S ≤ jacobson ⊥`）を足して'
  + '第 1366 `isDiscreteValuationRing_integralClosure`'
  + '（**完備 DVR の有限分離整閉包は DVR**）を得た。'
  + '★第 1367-1368 で下流が要求する形（`IsAdicComplete (maximalIdeal C) C`・`IsFractionRing C L`）に直し、'
  + '第 1369 で分岐指数 `e`（`m_A C = m_C^e`、`1 ≤ e ≤ [L:K]`）を取り出した。'
  + '★★第 1370-1372 で「不分岐の仮定 `hp`」を `e` 倍版 `hpe` に一般化した'
  + '——`exists_h2_h1_of_bad_prime` の証明で `hp` が使われるのは'
  + '`vAdd_algebraMap_eq_valAdd` の**ただ 1 箇所**だったので、そこだけ直せば道が通った'
  + '（`exists_h2_h1_of_bad_prime_ram`、第 1372）。';
e.remaining2026_09_02 = '★残るのは**組み立ての配管**であって、mathlib の欠落ではない:'
  + '(1) `L_v′/L_v` の付値の延長公式 `v_{L_v′}(y) = v_{L_v}(y)^e`（`y ∈ L_v`）——'
  + '第 1369 の `m_A C = m_C^e` から出るはず。'
  + '(2) 拡大の上での `IsMinimal`・`IsElliptic`・`HasSplitMultiplicativeReduction`——'
  + '不分岐版は `Found/GaloisRep/UnramQuad.lean`（第 1029-1033）にあり、'
  + '分岐版は付値の順序が保たれるので同様に出る。'
  + '(3) `ζ_l ∈ L_v′` の構成（`L_v(ζ_l)` は `≤ l−1` 次なので `l ∤ e`）。'
  + '(4) `SSCurve` の言葉（`E.alg ≃ₐ[E.fld] E.alg`）への翻訳。';
fs.writeFileSync(p, JSON.stringify(j, null, 2) + '\n', 'utf8');
console.log('resolved + remaining recorded');
