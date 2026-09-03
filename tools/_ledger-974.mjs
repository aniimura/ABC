// 2026-09-01（第 974）—— Vélu の商は変数変換と可換。界面の測定も記録。
// ★CLAUDE.md の規約により、台帳の更新は node -e ではなく .mjs で書く。
import { readFileSync, writeFileSync } from 'node:fs';

const P = 'D:/Math_ABC3/ResearchPaper/mathlib-gap.json';
const d = JSON.parse(readFileSync(P, 'utf8'));

d.lemma35VeluVcTransport20260901 = {
  measuredAt: '2026-09-01',
  block: '第 974',
  supersedes: 'lemma35DeltaFromJExp20260901',
  metric: 'GenEll 必要分 18 / 24（§1 9/9・§2 0/1・§3 6/9・§4 3/5）',
  closedThisStretch: [
    '★★★★★★★★★★★★`veluQuotientFull_vcPoint_eq`（974）——' +
      '`E′ = E/⟨Q⟩` なら `C • E′ = (C • E)/⟨C • Q⟩`。' +
      '☆点集合を 949 で捻り、`veluQuotientFull_variableChange` を当てるだけ。',
  ],
  interfaceFinding:
    '★★★★**界面の測定（重要）**: mathlib の `IsMinimal R W` は' +
    '「恒等変換が Δ の付値を最大にする」——すなわち **`W` 自身が極小モデル**という' +
    '主張であって、変数変換で不変ではない。' +
    '☆一方 `SemistableAt p E` が与えるのは「ある `C` の後で極小」である。' +
    '★したがって第 972 に渡せるのは `E`・`E′` そのものではなく `C • E`・`C′ • E′` である。' +
    '☆`minDeltaExp` は変数変換で不変（`minDeltaExp_variableChange`、在庫）なので' +
    '結論は戻せるが、**`C` と `C′` が別々**なので Vélu の関係 `hE′` が直接は渡らない。' +
    '★解決の道は 2 つ: ' +
    '(α) 974 で `C • E′` に移し、`C • E′` と `C′ • E′` の差は `j` が同じことで吸収する' +
    '（第 967 の `hE′` を `j` の等式に緩める）、' +
    '(β) 第 892/972 を「極小モデルを別データとして受ける」形に緩める。',
  remaining:
    '☆残るのは上の (α) か (β) を選んで界面を整えたうえで、' +
    '各悪い素点の局所データ（908・909・963・964・953・932・961・962・973）を並べる段。',
  note:
    '★944-974 の 31 ブロック。`Theorem 2.1`（§2）だけは独立で [NCBelyi] の形式化を要する。',
};

writeFileSync(P, JSON.stringify(d, null, 2) + '\n', 'utf8');
console.log('mathlib-gap.json: lemma35VeluVcTransport20260901 を書いた');
