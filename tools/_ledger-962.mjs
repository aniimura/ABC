// 2026-09-01（第 962）—— (d5)：Vélu の商の Δ は l¹² 倍、したがって楕円。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35VeluDelta20260901 = {
  measuredAt: '2026-09-01',
  block: '第 962',
  supersedes: 'lemma35VeluWMu20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★`Delta_velu_tate_eq`（962）——`Δ(Vélu の商) = l¹²·Δ(E_{q^l})`。' +
      '☆`1728Δ = c₄³ − c₆²` と `c₄ ↦ l⁴c₄`（第 884）・`c₆ ↦ l⁶c₆`（第 885）だけ。' +
      '`R` は `CharZero` の整域なので `1728` は消せる。',
    '★★★★★★★★★★★★★★★★`isElliptic_veluCurve_tate_map`（962）——' +
      '`l` が単元で `E_{q^l}` が楕円なら Vélu の商も楕円。' +
      '☆**同種写像の理論は一切要らない**——`c₄`・`c₆` の等式だけで済む。',
  ],
  remainingLayers: [
    '☆(D3) の残り: (a) `E ⊗ Lv` の分裂性、(b1) `hp`（完備化の付値の橋）。' +
      '★(d) は 956-962 の 7 ブロックで完結した。',
  ],
  note:
    '★944-962 の 19 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35VeluDelta20260901 を書いた');
