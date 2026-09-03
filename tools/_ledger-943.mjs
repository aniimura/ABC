// 2026-09-01（第 941-943）—— Lemma 3.5 に残る 1 節点と、その中の「拡大の層」。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35Node20260901 = {
  measuredAt: '2026-09-01',
  block: '第 941-943',
  supersedes: 'lemma35Assembly20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  singleNode:
    '★`Skeleton/GenEll/TateLocalModel.lean` の `isMuAtBadPrimes_of_veluQuotient`（941）が' +
    '`Lemma 3.5` に残るただ 1 つの節点である。' +
    '`lemma_3_5_velu_bad_delta`(903, 証明済み) がその結論 `IsMuAtBadPrimes` を' +
    'そのまま消費する。',
  correction:
    '★★第 941 で「残るのは機械的な組み立てだけ」と書いたのは不正確だった（第 942 で訂正）。' +
    '`lemma_3_2_i_tate_all`（906 が使う）は `K ⊆ L`（`IsGalois K L`）を受けている——' +
    '`μ_l` や `l`-捩れ点は完備化 `Lv` の中にあるとは限らない。' +
    '☆したがって **`Lv` の有限 Galois 拡大の層**が別に要る。',
  newTool:
    '★`isAdicComplete_valuationSubring`（943）——第 897 の一般形。' +
    '任意の完備な離散付値体の整数環で `IsAdicComplete` が取れる。' +
    '☆拡大の層で必要になるので先に用意した。',
  remainingLayers: [
    '☆(A) `Lv` の有限 Galois 拡大 `L′`（`μ_l` と `l`-捩れを含む）とその整数環を立てる。' +
      '`IsAdicComplete` は 943 が与える。',
    '☆(B) `Δ_min` が拡大でどう変わるか（分岐指数）を扱う層。',
    '☆(C) `H` の点が `L′`-有理になったうえで、Vélu の商の係数が `Lv` に降りること' +
      '（Vélu の式は対称式）。★本プロジェクトの `lemma_3_5_velu` 系は ' +
      '`Q : E.toAffine.Point`（`L`-有理）で受けているので、その逸脱も記録済み。',
    '☆(D) 上を並べて `isMuAtBadPrimes_of_veluQuotient` を閉じる。' +
      '分裂/非分裂の連鎖はすでに全部ある（906/927/904、938/925/940/926/929）。',
  ],
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35Node20260901 を書いた');
