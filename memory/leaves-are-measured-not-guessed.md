---
name: leaves-are-measured-not-guessed
description: 葉はページ順では分からない。node tools/frdi-leaves.mjs で「未実装依存の数」と「波及」を測ってから着手する。
metadata:
  type: feedback
---

進め方は「葉から形式化する」だが、★★**どれが葉かは物理ページ順では分からない**。

`node tools/frdi-leaves.mjs` は、残っている項目それぞれについて
「原文が番号で挙げた依存のうち**まだ実装していない**ものの数」を数え、

1. **依存 0 の項目(＝葉)**
2. **その葉を開けたときの波及**（何件の未実装依存が減るか）

を出す。

**Why:** 2026-08-18 の実測で、残り 25 件のうち葉は 10 件。波及は
`Theorem 5.2` が 6 件、`Proposition 4.4` が 4 件だった。
★ページ順に見ていた限り `Theorem 5.2` は「§5 の大物」で後回しに見えたが、
実際には**依存 0 の葉**で、`Proposition 5.3`・`Example 6.1`・`Example 6.3`・
`Theorem 6.2`・`Proposition 5.5`・`Theorem 6.4` がそこで詰まっていた。

**How to apply:**

1. 着手対象を選ぶ前に `node tools/frdi-leaves.mjs` を回す。
2. 葉のうち**波及の大きいもの**から取る。波及 0 の葉（`Remark 4.9.1` 等）は、
   カウンタを 1 動かす以上の意味はない——急ぐときだけ。
3. ★**限界**: pdftotext 経由なので番号の無い依存（§0 の語彙、「immediate」の段）は
   数えられない。当たりを付けるためだけに使い、実装で出た必要物は
   [[leaf-first-with-graph-feedback]] のとおりスケルトンに足してグラフを更新する。

関連: [[leaf-first-with-graph-feedback]] [[lean-check-via-mcp-repl]]
