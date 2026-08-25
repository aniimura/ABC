---
name: stale-status-read-lean-first
description: 「残っている」と台帳が書いた条は、着手前に必ず Lean の側を読む。6 回とも既に閉じていた。
metadata:
  type: feedback
---

`ResearchPaper/frdi-decomposition.json` の `status` / `note` は**書いた時点の見立て**であり、
その後に閉じても更新されないことがある。★**着手前に必ず Lean の側を読むこと。**

実測の手順(どれも 1 分以内):

1. `node tools/decl-index.mjs` で `.cache/decl-index.txt` を作り、節点の `decl` を grep
2. 対応するファイルの `grep -c sorry`
3. `mcp__abc3-lean__lean_check` に `#check @Foo` を並べて、在庫の宣言名を確定

**Why:** 同じ型の見立て違いが **6 回**起きた:

| 回 | 節点 | 実際 |
|---|---|---|
| 1 | `Prop 5.3` の右の縦の矢印 | 済んでいた |
| 2 | `hrefl` | 済んでいた |
| 3 | `hprim` | 済んでいた |
| 4 | `p55iii-pf` | 済んでいた |
| 5 | `thm62-iii-ratstd` | 「残るのは配線」と書いてあったが配線も済んでいた |
| 6 | `ex63-model`(2026-08-25) | ★`Skeleton/Divisor/ArithDivisor.lean` の `sorry` 11 個は**実装済み `Found` の古い写し**。算術 Frobenioid は `Found/Divisor/` の 17 ファイル・4,291 行に `sorry` 0 で組み上がっていた |

★★6 回目は特に危険な形だった —— **`Skeleton` に `sorry` があるのに `Found` では閉じている**。
`grep sorry` だけで判断すると「未着手」と読んでしまう。
★`Skeleton` は「壁の見取り図」なので、閉じた後も委譲に書き換えられていないことがある。

**How to apply:**

- 節点に着手する前に `decl` の grep と `#check` を必ず走らせる。
- `Skeleton` に `sorry` があっても、同名・同趣旨の宣言が `Found` に無いか探す。
- 閉じたら `Skeleton` を**委譲に書き換える**(`:= ABC3.Found.…`)。そうすれば次に
  `grep sorry` した人が誤読しない。★このセッションで `Thm62Slim` / `Thm64Deg` /
  `Divisor/Normalization` / `Divisor/Hartogs` はそう直した。
- 台帳の `note` は「測った日付つき」で追記する。古い段落は消さず、
  「★★★★訂正(日付)」の見出しで**上書きせずに並べる**(判断の履歴が残る)。

関連: [[leaves-are-measured-not-guessed]] / [[measure-mathlib-before-skeleton]] /
[[frdi-split-nonisotropic-not-derivable]]
