// 2026-09-01（第 968）—— 曲線の等式に沿って点を運ぶ。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35PointImageEq20260901 = {
  measuredAt: '2026-09-01',
  block: '第 968',
  supersedes: 'lemma35JMapVelu20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★★★★★`exists_point_image_eq`（968）——' +
      '曲線が等しければ、位数 `l` の点とその**座標集合**を運べる。`subst` 一発。',
  ],
  designNote:
    '★★点の輸送を `▸` で書き、`veluQuotientFull` の楕円性インスタンスまで' +
    '`∧` の中で運ぼうとすると詰まる——`∧` の第 2 項は第 1 項をインスタンスとして使えない。' +
    '☆そこで**インスタンスを含まない形**（座標集合の一致）で出した。' +
    '呼ぶ側はそれで `rw` すれば楕円性も `j` の等式も一緒に移る。',
  assemblyPlan: [
    '☆`isMuAtBadPrimes_of_veluQuotient` の本体の道筋（部品はすべて在庫）:',
    '(1) `p` を固定し `jExp p E < 0` を仮定。`Lv := p.adicCompletion L`、`R := 整数環`。',
    '(2) 954 で `C`・`hC`・`hc4`（`E` と `E′` の両方）。',
    '(3) 909 で乗法還元、963 で分裂性（非分裂なら 938 で捻り、最後に 929 で降ろす）。',
    '(4) 964 で `hp`。953 で `hlu`。',
    '(5) `tateParamR_spec` ＋ 944（`tateModel_baseChange`）で ' +
      '`(E_q) ⊗ Lv = C_K • (E ⊗ Lv)`。',
    '(6) 967 で `hW′j` を `C_K • (E ⊗ Lv)` の形で得る。',
    '(7) 968 でそれを Tate モデルの形に移し、点 `P` を得る。',
    '(8) 961 ＋ 962 で `hvw`。932 で `hcop`。',
    '(9) 965 に渡して `minDeltaExp p E′ = l · minDeltaExp p E`。',
  ],
  note:
    '★944-968 の 25 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35PointImageEq20260901 を書いた');
