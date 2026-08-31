// 2026-09-01（第 966）—— 捉れ点を Tate モデルまで運ぶ。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35PointTransport20260901 = {
  measuredAt: '2026-09-01',
  block: '第 966',
  supersedes: 'lemma35LocalTerminal20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★`nsmul_eq_zero_of_addOrderOf` / `ne_zero_of_addOrderOf_prime`（966）——' +
      '`addOrderOf Q = l` から `l • Q = 0` と `Q ≠ 0`。',
    '★★★★`rhPoint_nsmul_eq_zero` / `rhPoint_ne_zero`（966）——体拡大で運ぶ。',
    '★★★★`vcPoint_nsmul_eq_zero` / `vcPoint_ne_zero`（966）——変数変換で運ぶ。',
    '☆これで `L` の上の `Q`（`addOrderOf Q = l`）から、' +
      'Tate モデルの上の `P`（`l • P = 0`、`P ≠ 0`）が作れる——' +
      '第 965 が受ける形である。' +
      '★位数がちょうど保たれることまでは要らないのが軽い。',
  ],
  plumbing:
    '★配管の実測: `rhPoint_nsmul` などの `•` は宣言側が `open scoped Classical` なので' +
    'Classical の `DecidableEq` から来る加法群の `SMul` を使っている。' +
    '呼ぶ側で `[DecidableEq F]` を束縛すると別インスタンスになり `rw` が当たらない。' +
    '☆`tools/lean-idioms.md` に追記した。',
  remainingLayers: [
    '☆残るのは `isMuAtBadPrimes_of_veluQuotient` の本体で、各悪い素点について' +
      '965 の仮説を並べる段。部品はすべて揃っている' +
      '（963/909/929・964・954・932・953・961・962・950・966）。',
  ],
  note:
    '★944-966 の 23 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35PointTransport20260901 を書いた');
