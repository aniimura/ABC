// 2026-09-01（第 972）—— 大域の Vélu の商で受ける形。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35GlobalVelu20260901 = {
  measuredAt: '2026-09-01',
  block: '第 972',
  supersedes: 'lemma35TateSideElliptic20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★★★★★`minDeltaExp_eq_mul_of_globalVelu`（972）——' +
      '第 965 は Tate モデルの上の点 `P`・`j` の一致 `hW′j`・商の楕円性 `hellQ` を受けていた。' +
      '☆本定理はその 4 つを**大域のデータ `Q`・`hQ`・`hE′` から作る**:' +
      ' `P` 系は 970、`ζ`・`uζ` は 947、`hellQ` は 971。',
  ],
  remaining:
    '★残るのは `isMuAtBadPrimes_of_veluQuotient` の本体で、' +
    '**各悪い素点で 972 の局所データを供給する段**だけである:' +
    ' `h`・`h′`（963＋909、非分裂は 938→929）、`hp`（964）、' +
    ' `C`・`hC`・`hc4`（954）、`hΔ`・`hcop`（932）、`hlu`（953）、' +
    ' `hql`・`h2`・`h2K`、`hvw`（961＋962）、' +
    ' そして `E ⊗ Lv`・`E′ ⊗ Lv` の極小性（908）。',
  note:
    '★944-972 の 29 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35GlobalVelu20260901 を書いた');
