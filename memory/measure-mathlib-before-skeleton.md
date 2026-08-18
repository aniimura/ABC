---
name: measure-mathlib-before-skeleton
description: 新しい債務に着手する前に、引用先が mathlib/FLT にあるか実測してから欠落だけをスケルトン化する
metadata:
  type: feedback
---

新しい obligation に着手するときは、**引用先を列挙 → mathlib/FLT の在庫を実測(grep + ファイルを開いて TODO を読む)→ 欠落だけをスケルトン化して依存グラフに節点追加 → 葉から形式化** の順を守る。

**Why:** B1(2026-08-18 達成、146 ブロック)で見積もりが膨らんだ原因は、引用先(Hartshorne II.5 相当)が mathlib にあると**仮定した**ことだった。実際には 2 つが TODO のまま残っており、自前で 27 ブロック(B1 全体の 19%)を要した:

- `tilde` がテンソル積を保つ(mathlib は「同値な命題は未解決」と TODO に明記)→ 17 ブロック
- 準連接 ⟹ `fromTildeΓ` 同型(`Mathlib/AlgebraicGeometry/Modules/Tilde.lean` 547 行の TODO)→ 10 ブロック

どちらも**ファイルを開けば 1 分で分かった**。対照的に、§9-149 で在庫を実測してから見積もった QC 群は、見積もり 26 に対し実績 12 で収まった。

**How to apply:** 単位は「引用した教科書」ではなく「**mathlib に無いもの**」。教科書全体をスケルトン化すると過大計上になる(大半は mathlib にある)。grep で名前を探すだけでなく、**見つかったファイルの TODO コメントを読む**こと——宣言があっても中身が未完のことがある。[[verify-insertion-not-just-ok]] と同じ「ok を信じず中身を見る」姿勢。
