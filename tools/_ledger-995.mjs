// 2026-09-01（第 995）—— E′ の側は j だけでよい（界面の大きな軽量化）。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35EprimeJOnly20260901 = {
  measuredAt: '2026-09-01',
  block: '第 995',
  supersedes: 'lemma35NormalizeU1_20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★★★★★`jExp_eq_mul_of_tateParam_pow`（995）——' +
      '`E′` の `j` が `E_{q^l}` の `j` なら `jExp p E′ = l · jExp p E`。' +
      '☆第 932 を `E` と `E′` の両方に当て、`vAdd (q^l) = l · vAdd q` で割るだけ。',
    '★★★★★★★★★★★★`minDeltaExp_eq_mul_of_jExp_mul`（995）——' +
      '`jExp` が `l` 倍なら `Δ_min` も `l` 倍（半安定・悪い素点で）。',
  ],
  bigSimplification:
    '★★★★**界面が大きく軽くなった**: 第 892／972 は `E′` についても' +
    '分裂乗法還元・極小モデル・`C′`・`hc4′` を要求していた。' +
    '☆だが `minDeltaExp = max(0, −jExp)`（在庫）なので、' +
    '**`jExp p E′ = l · jExp p E` さえ出れば済む**。' +
    'そして `jExp` は第 932 により **`j` の一致だけ**から出る。' +
    '★したがって `E′` の側に分裂性も極小モデルも要らない——' +
    '要るのは「`E′ ⊗ Lv` の `j` が `E_{q^l} ⊗ Lv` の `j` に等しい」ことだけである。',
  remaining:
    '☆残るのは (i) `E` の側の局所データを揃える段（954・976・993・994、すべて在庫）、' +
    '(ii) `j` の一致を 950／967／970 から取り出す段、' +
    '(iii) 995 に流して `IsMuAtBadPrimes` を閉じる段。',
  note:
    '★944-995 の 47 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35EprimeJOnly20260901 を書いた');
