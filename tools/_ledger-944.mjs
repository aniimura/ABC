// 2026-09-01（第 944）—— 分岐指数の層（第 943 の (B)）が消えた。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35Descend20260901 = {
  measuredAt: '2026-09-01',
  block: '第 944',
  supersedes: 'lemma35Node20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★`tateCurveAt_map_comm`（944）——`E_{σq} ⊗ K′ = (E_q ⊗ K) ⊗ K′`。' +
      '`tateCurveAt_map`（869）と `map_map` だけ。',
    '★`tateModel_baseChange` / `j_tateModel_eq`（944）——' +
      'Tate モデルを `K` に戻すと変数変換を除いてもとの曲線、したがって `j` も同じ。',
    '★★★★★★★★★★★★★★★★`tateParamR_map`（944）——**`q_{W ⊗ K′} = σ(q_W)`**。' +
      'Tate 母数は体の拡大で変わらない。',
    '★★★★★★★★`tateParam_descend`（944）——拡大の上で得た `q′ = q^l` を' +
      '`σ` の単射性でそのまま `R` に降ろす。',
  ],
  layerBRemoved:
    '★★★★第 943 で「残る層」に挙げた **(B) `Δ_min` が拡大でどう変わるか（分岐指数）は不要**になった。' +
    '☆見方を変える: 拡大の上で得た結論を `Δ_min` の段で降ろすのではなく、' +
    '**Tate 母数の段で降ろす**。`q` は `R` の元であり `σ : R → R′` は単射だから、' +
    '`σ(q_{E′}) = σ(q_E)^l` からその場で `q_{E′} = q_E^l` が出る。' +
    '`Δ_min` はその後で `R` の中だけで計れるので、分岐指数は一切現れない。',
  remainingLayers: [
    '☆(A) `Lv` の有限拡大 `L′`（`μ_l` と `l`-捩れを含む）とその整数環を立てる。' +
      '`IsAdicComplete` は 943 が与える。' +
      '★`tateParamR_map` が要求するのは `σ`（局所準同型）・`φ`（体準同型）・可換性だけ。',
    '☆(C) `H` の点が `L′`-有理になったうえで、Vélu の商の係数が `Lv` に降りること' +
      '（Vélu の式は対称式）。',
    '☆(D) 上を並べて `isMuAtBadPrimes_of_veluQuotient` を閉じる。' +
      '分裂/非分裂の連鎖はすでに全部ある（906/927/904、938/925/940/926/929）。',
  ],
  note:
    '★外部文献への未証明の依存は無い。ht^Falt は第 704 で内製、' +
    'Néron–Ogg–Shafarevich は第 902 で不要、' +
    '形式逆関数定理は第 875（単射）＋第 933（全射）で両側が揃った。' +
    '★`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35Descend20260901 を書いた');
