// 2026-09-01（第 946）—— 有理点なら Galois は不要。組み立ての道具が揃った。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35Rational20260901 = {
  measuredAt: '2026-09-01',
  block: '第 946',
  supersedes: 'lemma35MuDescend20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★`dvd_of_pow_eq_zpow_of_coprime`（946）——`x^l = Q^k` と `l ∤ v(Q)` から `l ∣ k`。' +
      '`l·v(x) = k·v(Q)` の 1 行。',
    '★★★★★★★★`exists_mu_of_rational`（946）——`K`-有理な位数 `l` の類は `μ_l` の類。',
    '★★★★★★★★★★★★`exists_mu_point_of_rational`（946）——点の側で述べた形。',
    '★★★★★★★★★★★★★★★★★★★★`exists_mu_point_dvr`（946）——' +
      '**完備な DVR だけで述べた形**。引数は `q`・`Δ ≠ 0`・`l` が素・`l ∤ v(q)`・`l • P = 0` の 5 つだけ。',
  ],
  galoisRemoved:
    '★★★★★★★★**`Lemma 3.2, (i)` が Galois 安定性を使うのは、点が拡大体にしかない場合に' +
    '`x` を握るためである**。`Lemma 3.5` の `Q` は `L`-有理だから `x ∈ Lvˣ` そのものであり、' +
    '`x^l = q^k ⇒ l·v(x) = k·v(q)` という**付値の計算 1 行**で `l ∣ k` が出る。' +
    '☆これで組み立てから `IsGalois`・拡大体 `L`・`hσv`・`q₀`・`v` が**すべて消えた**。' +
    '★Tate 一意化 `Φ` は `dvrTatePhiAddEquiv`（899）が無条件に与える。',
  remainingLayers: [
    '☆(D) 組み立てのみ。残るのは 3 点:' +
      ' (D1) `hcop : True` を本物（`l ∤ v_p(q_E)`、原文の「l is prime to the local heights」）にする、' +
      ' (D2) `H = ⟨Q⟩` の点集合が `{[ζ^i]}` であることを `pointCoords` の水準で言う、' +
      ' (D3) 927 → 904 → 903 と流す。' +
      '★体論的な障害はもう無い。',
  ],
  note:
    '★第 943 の (A)(B)(C) は 945・944・945 で消え、Galois そのものが 946 で消えた。' +
    '外部文献への未証明の依存は無い。' +
    '★`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35Rational20260901 を書いた');
