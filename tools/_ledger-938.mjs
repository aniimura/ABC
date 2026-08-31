// 2026-09-01（第 935-938）—— 非分裂の降下の (1) が閉じた。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35Split20260901 = {
  measuredAt: '2026-09-01',
  block: '第 935-938',
  supersedes: 'lemma35Inverse20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★`integralModel_eq_of_map_eq` / `integralModel_quadTwist`（935）——' +
      '整モデルは一意（底変換が一致する R 上の曲線は係数ごとに一意）。' +
      'したがって「捻りの整モデルは整モデルの捻り」。',
    '★`exists_sq_of_twist_of_not_isSquare`（936）——' +
      '非分裂の 2 次式は非平方で捻ると根をもつ。',
    '★`quadTwist_splitConst` / `quadTwist_splitLin`（937）——' +
      'Splits の 2 次式は c₄ ↦ d²c₄、定数項 ↦ d³×、X の係数は 0。',
    '★★★★`splits_quadTwist_of_not_isSquare`（938）——' +
      'mathlib の `HasSplitMultiplicativeReduction` の `Splits` フィールドそのものを' +
      '捻った整モデルについて出す。**残る 3 段のうち (1) はこれで閉じた**。',
  ],
  remaining: [
    '☆(2) 捻った対の `hv`・`hw`（μ_l の和との同定）——分裂の場合と同じ帳簿。' +
      '925 が曲線の等式（捻った商は veluCurve の形）を、926 が変換則を与えている。',
    '☆(3) 組み立て: 938 で捻りが分裂と分かるので、捻った対に分裂の連鎖' +
      '（906 → 927 → 904）を当て、929 で降ろし、903 に流す。',
  ],
  note:
    '★外部文献への未証明の依存は無い。ht^Falt は第 704 で内製、' +
    'Néron–Ogg–Shafarevich は第 902 で不要になり、' +
    '形式逆関数定理は第 875（単射）と第 933（全射）で両側が揃った。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35Split20260901 を書いた');
