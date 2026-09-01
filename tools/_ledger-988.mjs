// 2026-09-01（第 988）—— 剰余体の有限性の核：𝓞 L は p を法として代表を尽くす。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35IntegerCongr20260901 = {
  measuredAt: '2026-09-01',
  block: '第 988',
  supersedes: 'residueFieldFiniteGap20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★`exists_integer_congr`（988）——' +
      '`v_p(y) ≤ 1` なる `y ∈ L` に対し `a ∈ 𝓞 L` を `v_p(y − a) < 1` に取れる。' +
      '☆第 985・987 で立てた道筋 (a)-(d) のうち **(b)(c) が閉じた**——' +
      'すなわち剰余体の有限性の**核**である。',
  ],
  howItWasDone: [
    '☆道は 4 段:',
    '(1) `y` の分母に現れる素点は有限個（`HeightOneSpectrum.Support.finite`）',
    '(2) `v_p(y) ≤ 1` だからそれらは `p` と異なる',
    '(3) 中国剰余（`IsDedekindDomain.exists_forall_sub_mem_ideal`）で ' +
      '「`p` では `1`、分母の素点では十分高い冪」に合同な `s` を取る',
    '(4) `s·y` は全付値が `≤ 1` なので整（`mem_integers_of_valuation_le_one`）、' +
      '`y − s·y = (1 − s)·y` の `p`-付値は `< 1`',
  ],
  remaining:
    '☆残るのは (i) `L` の稠密性（`denseRange_algebraMap`、mathlib にあり）と 988 を繋いで ' +
    '`𝓞 L → ResidueField (p.adicCompletionIntegers L)` の全射性を出し、' +
    '`Ideal.finiteQuotientOfFreeOfNeBot` で `Finite` を得る段、' +
    '(ii) `C • E` の正規化、(iii) 捻り `d` の降下、(iv) 972 に並べる段。',
  note:
    '★944-988 の 41 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35IntegerCongr20260901 を書いた');
