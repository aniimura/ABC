---
name: mathlib-index-nonascii-truncation
description: decl-index の名前欄は非 ASCII で切れていた(preΨ₄ → pre)。mathlib 24,659 宣言が誤った名前で載っていた。2026-09-06 に修正
metadata:
  type: reference
---

**`tools/decl-index.mjs` の `NAME` が ASCII の `\w` だけで組まれていたため、
索引の名前欄が最初の非 ASCII 文字で切れていた。**

実測（2026-09-06、修正前）:

- `WeierstrassCurve.preΨ₄` → 索引上 `WeierstrassCurve.pre`
- `coeff_Ψ` → 索引上 `coeff_`
- `Ψ` で**名前欄**を grep すると **0 件**。しかし statement 欄には 133 行あった

修正後の測定:

| | 修正前 | 修正後 |
|---|---|---|
| mathlib 索引の行数 | 247,821 | **249,273**（+1,452 は行ごと落ちていた） |
| ABC3 索引の行数 | 23,107 | 23,170 |
| 名前欄が `Ψ` を含む行 | **0** | **111** |
| 名前欄が非 ASCII を含む行(mathlib) | — | **24,659**（≒1 割が誤った名前で載っていた） |

**Why:** これは [[mathlib-frobenius-element-exists]] を含む「不在の誤判定」5 件と
**同じ回路**である（`ULift.field` / `continuousCohomology` / `Ẑ` / `CompactSpace Gal` /
`IsArithFrobAt`）。無名 instance を索引に入れた 2026-09-05 の拡張と動機も同じ。
★見つけたのは Weil 対の在庫調査で、`WeierstrassCurve.Ψ/Φ` を「無い」と誤判定しかけた瞬間である。

**How to apply:**

- ★**索引の名前欄が `0 件` を返しても、非 ASCII を含む名前なら疑う**。
  statement 欄（4 列目）は切れていないので、そちらでも引く。
- 入れた文字類は Latin Extended / ギリシャ / 音声拡張 / 上下付き / 文字風記号 / 代理対。
  ★**矢印（`←-⇿`）と論理記号は入れていない** —— 入れると
  `def f : A → B` の型まで名前に食い込む。回帰試験の考え方はこの 2 例で足りる。
- 「無い」と書く前に `node tools/absent-recheck.mjs --try '<正規表現>'`。

★あわせて見つけた `check.mjs` の脆さ:
**G6 の obligation 区切りは文字列の中まで走査する**ので、docstring や `.needs` の
本文に `.implicitStep` のような**先頭ドット付きの綴り**を書くと区切りが誤爆し、
`stripStr` が壊れて本文中の数値が頁番号として拾われる（実際に「物理 p.1036 は範囲外」
という偽の NG が出た）。→ 記録本文では先頭のドットを外して書く。
