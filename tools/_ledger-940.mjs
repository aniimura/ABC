// 2026-09-01（第 939-940）—— 非分裂の降下は残り (3) の組み立てだけになった。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35Assembly20260901 = {
  measuredAt: '2026-09-01',
  block: '第 939-940',
  supersedes: 'lemma35Split20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★`veluV2_quadTwist` / `veluU_quadTwist`（939）——`a₁ = 0` の形では ' +
      '`v_Q = 3x² + 2a₂x + a₄` で **`y` に依らない**。捻りで `d²` 倍。',
    '★`veluVFull_quadTwist`（940）——捻った点集合の `v` の和は `d²` 倍。' +
      '`y` に依らないので捻った点の `y` 座標（`√d` を含む）を一切見ずに済む。',
    '☆これで残る 3 段のうち **(1) と (2) が閉じた**。',
  ],
  remaining: [
    '☆**(3) の組み立てだけ**である。並べる順序:' +
      ' 938（`d` を非平方単数に取れば捻りは `p` で分裂乗法還元）' +
      ' → 935（捻りの整モデルは整モデルの捻り）' +
      ' → 925（捻った商は `veluCurve` の形）＋ 940（`v` の和）' +
      ' → 926（Tate モデルへ移す）' +
      ' → 906 / 927 / 904（分裂の連鎖を捻った対に当てる）' +
      ' → 929（もとの対へ降ろす）' +
      ' → 903（`Lemma 3.5` へ流す）。',
  ],
  note:
    '★外部文献への未証明の依存は無い。ht^Falt は第 704 で内製、' +
    'Néron–Ogg–Shafarevich は第 902 で不要、' +
    '形式逆関数定理は第 875（単射）＋第 933（全射）で両側が揃った。' +
    '☆`Lemma 3.5` が閉じれば `Lemma 3.7` → `Theorem 3.8` → ' +
    '`Corollary 4.3` / `4.4` は既存の導出（第 367）が動く。' +
    '★`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35Assembly20260901 を書いた');
