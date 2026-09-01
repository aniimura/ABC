// 2026-09-01（第 990）—— 剰余の代表は 𝓞 L から取れる（捻り d の降下）。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35ResidueRep20260901 = {
  measuredAt: '2026-09-01',
  block: '第 990',
  supersedes: 'lemma35ResidueFieldFinite20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★`exists_integer_residue_eq`（990）——' +
      '完備化の整数環の剰余の代表は `𝓞 L` から取れる（第 989 の全射性の中身）。' +
      '☆第 982 が出す捻り `d` は `R` の元だが、第 929 が受けるのは `L` の元である。' +
      '★2 次式の係数は `d` に `φ d` を通してしか依存しないので、' +
      '剰余が同じ `𝓞 L` の元に取り替えれば `Splits` はそのまま移る。',
  ],
  remaining:
    '☆残るのは (i) `C • E` を `a₁ = a₃ = 0` に正規化して `p` での整性を保つ段、' +
    '(ii) 982 の二者択一を `hasSplitMultiplicativeReduction_baseChange` に流す段' +
    '（990 で `d` を `L` に降ろしてから 929 へ）、' +
    '(iii) 972 に並べて `minDeltaExp_variableChange` で戻す段。' +
    '★mathlib の穴はもう無い——すべてプロジェクト内の配管である。',
  note:
    '★944-990 の 43 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35ResidueRep20260901 を書いた');
