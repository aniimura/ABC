// 2026-09-01（第 971）—— Tate モデル側の商の楕円性。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35TateSideElliptic20260901 = {
  measuredAt: '2026-09-01',
  block: '第 971',
  supersedes: 'lemma35TateModelPoint20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★`isElliptic_veluQuotient_tate_mu`（971）——' +
      '`μ_l` による Vélu の商（**Tate モデル側**）は楕円。' +
      '☆出どころは `hvw` が持っている `veluCurve` の楕円性で、' +
      '`veluQuotientFull_tate_mu`（890）が両者を繋ぐ。' +
      '★これで第 965 の `hellQ` が `hvw` から直に出る。',
  ],
  twoSides:
    '★楕円性の producer が 2 つ揃った: ' +
    '969 は `C • (E ⊗ Lv)` **側**の商について（`E′` から来る）、' +
    '971 は **Tate モデル側**の商について（`hvw` から来る）。' +
    '☆二つは 968 で移り合う。',
  remainingLayers: [
    '☆残るのは 970 ＋ 947 ＋ 971 ＋ 965 を 1 つの証明に並べる段と、' +
      '各悪い素点で局所データ（954・909・963・929・964・953・932・961・962）を供給する段。',
  ],
  note:
    '★944-971 の 28 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35TateSideElliptic20260901 を書いた');
