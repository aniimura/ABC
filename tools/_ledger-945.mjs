// 2026-09-01（第 945）—— 拡大の層（第 943 の (A)）も (C) も消えた。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35MuDescend20260901 = {
  measuredAt: '2026-09-01',
  block: '第 945',
  supersedes: 'lemma35Descend20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★`mu_fixed_of_pointMap_eq`（945）——' +
      '`tatePhi([ζ])` が σ で不変なら `σζ = ζ`。' +
      '`[σζ] = [ζ]` から `σζ = ζ·Q^k` を取り、付値で `k·v(Q) = 0`、`v(Q) > 0` で `k = 0`。' +
      '☆`ζ` が 1 の冪根であることさえ使わない。',
    '★★★★★★★★`exists_mu_base_of_rational`（945）——`ζ` は `K` から来る。' +
      'mathlib の `InfiniteGalois.mem_range_algebraMap_iff_fixed` を当てるだけ。',
    '★★★★★★★★★★★★★★★★`exists_mu_point_rational`（945）——第 906 と繋いだ形。' +
      '`K`-有理な位数 `l` の点は `K` の中の `μ_l` の点である。',
  ],
  layersACRemoved:
    '★★★★第 943 で「残る層」に挙げた **(A)（`μ_l` を含む有限 Galois 拡大 `L′` を立てる）は不要**になった。' +
    '☆`Lemma 3.5` の `Q` は `Q : E.toAffine.Point`、すなわち **`L`-有理な点**である。' +
    'したがって `tatePhi([ζ]) = Q ⊗ Lv` は Galois で不変であり、' +
    '上の降下で `ζ` は自動的に `Lv` に入る。' +
    '★★同じ理由で **(C)（Vélu の係数の降下）も消える**——' +
    '`H = ⟨Q⟩` の点ははじめから `L`-有理なので、Vélu の和は `Lv` の中で計算できる。',
  remainingLayers: [
    '☆(D) 組み立てのみ。`isMuAtBadPrimes_of_veluQuotient` の中で' +
      ' 945（`ζ ∈ Lv`）→ 927（`q_{E′} = q_E^l`）→ 904（`Δ_min` が `l` 倍）' +
      ' と並べ、非分裂は 938/925/940/926/929 で捻りに移してから同じ連鎖を当てる。' +
      '★体論的な障害はもう無い。',
  ],
  note:
    '★第 943 の (A)(B)(C) はそれぞれ 945・944・945 で消えた。' +
    '外部文献への未証明の依存は無い。' +
    '★`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35MuDescend20260901 を書いた');
