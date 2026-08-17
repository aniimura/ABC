---
name: frdi-three-walls
description: 「§2 以降を満点」は現在の在庫では到達不能。壁は 3 つで、node tools/frdi-blocked.mjs が節ごとの到達可能な上限を出す。
metadata:
  type: project
---

2026-08-18 の機械測定(`node tools/frdi-blocked.mjs`)。

★★**「§2 8/8・§4 10/10・§5 7/7・§6 5/5」は現在の在庫では到達できない。**

| 節 | 現在 | **到達可能な上限** |
|---|---|---|
| §1 | 13/15 | 15/15 |
| §2 | 7/8 | **7/8** |
| §3 | 7/9 | 9/9 |
| §4 | 3/10 | **5/10** |
| §5 | 0/7 | **2/7** |
| §6 | 0/5 | **3/5** |

**壁は 3 つだけ**(残り 24 件のうち 13 件を止めている):

1. **mathlib に pro-`l` 群が無い**(`Definition 2.8, (ii)`)—— §2 の最後の 1 件。
   位相的有限生成な副有限アーベル群の `M ≅ ∏_l M[l]`。原文 p.106–107 が等式として使う。
2. ★★**`Proposition 4.4, (ii)` の一般ケース**(`otriBase` ⟺ `𝒪^▷(A^birat)` の可換性)
   —— §4 で 5 件・§5 で 5 件。**最大の律速**。
   ★2026-08-18 に **model Frobenioid** と **birat-Frobenius-normalized 型**では閉じた
   (`Thm52Birat.lean` / `Prop44Otri.lean`)。残るのは一般の `𝒞`。
3. **mathlib に six exponentials theorem が無い**(`Lemma 6.5, (ii)`)—— §6 で 2 件。

**How to apply:**

- ★**壁が無い節は §1 と §3 だけ**。次に満点にできるのは **§3**
  (`Proposition 3.2` → `Theorem 3.4`、`Proposition 1.14` 経由)。
- 「§2 以降を満点」を目標に据え続けるなら、**壁 1・3 は mathlib への貢献**、
  **壁 2 は原典の「routine exercise」を実際に埋める研究**になる。どれを取るかは人が決めること。
- ★数字を確認するときは `node tools/frdi-blocked.mjs` と
  `node tools/frdi-leaves.mjs` を回す(手で数えない)。

関連: [[leaves-are-measured-not-guessed]] [[leaf-first-with-graph-feedback]]
