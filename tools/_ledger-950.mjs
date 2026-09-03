// 2026-09-01（第 949-950）—— (D2) の点集合の帳簿が閉じた。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35PointSet20260901 = {
  measuredAt: '2026-09-01',
  block: '第 949-950',
  supersedes: 'lemma35LocalStep20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★`pointCoords_neg` / `nsmul_eq_neg_nsmul_of_addOrderOf`（949）——' +
      '`-(k·Q) = (l-k)·Q`（`l·Q = 0` だから）。',
    '★★★★★★★★★★★★`pointCoords_image_negY_stable`（949）——' +
      '`⟨Q⟩∖{O}` の座標集合は反転で安定。' +
      '☆これが `j_veluQuotientFull_variableChange`（915）の `hS` そのものである。',
    '★★★★★★★★★★★★★★★★`image_pointCoords_vcPoint_nsmul`（949）——' +
      '`image_pointCoords_rhPoint_nsmul`（体準同型版）の**変数変換版**。',
    '★★★★★★★★★★★★★★★★★★★★`j_veluQuotientFull_nsmul_variableChange`（950）——' +
      '**(D2) そのもの**。`E ⊗ Lv` の側で与えられた Vélu の商を' +
      'Tate モデル `C • (E ⊗ Lv)` の側の商に移す。',
  ],
  remainingLayers: [
    '☆(D3) 各悪い素点で 948 → 904 → 903 と流す段だけ。' +
      '950 が `j` の一致を与え、948 が `q_{E′} = q_E^l` を与え、' +
      '904 が `Δ_min` の `l` 倍を与え、903 が `Lemma 3.5` に流す。' +
      '★非分裂は 938/925/940/926/929 で捻りに移してから同じ連鎖を当てる。',
  ],
  note:
    '★第 943 の (A)(B)(C) は 945・944・945 で消え、Galois は 946 で消え、' +
    '局所の段は 947-948 で 1 つの定理になり、点集合の帳簿は 949-950 で閉じた。' +
    '☆残るのは各悪い素点で既存の定理を並べる段だけである。' +
    '★`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35PointSet20260901 を書いた');
