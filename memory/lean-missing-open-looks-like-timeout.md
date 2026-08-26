---
name: lean-missing-open-looks-like-timeout
description: Lean 4 で `open` を書き忘れると「unknown identifier」ではなく whnf/isDefEq のタイムアウトとして現れることがある
metadata:
  type: feedback
---

`ABC3/Found/Divisor/NormNormal.lean` を書いたとき、`open ABC3.Found.FrdI` を
忘れたために `FinSub V.functionField Kbar` が解決できなかったが、
エラーは **`unknown identifier` ではなく
`(deterministic) timeout at whnf` / `at isDefEq`** として、
しかも**定理の本体ではなく型（statement）の途中の列位置**に出た。
`set_option maxHeartbeats` を上げても消えず、1 時間近く別の原因を疑って空回りした。

**Why:** 未解決の識別子があると、エラボレータは別の解釈（他の名前空間、
コアーション、instance 探索）を延々と試すため、
名前解決の失敗が**性能の症状**に化ける。

**How to apply:** Lean で「型の途中で whnf / isDefEq のタイムアウト」が出たら、
`maxHeartbeats` を上げる前に **まずその行に出てくる識別子がすべて解決しているか**を
確かめる（`#check` を 1 行足す、あるいは同じ識別子を使っている既存ファイルの
`open` 行と突き合わせる）。関連: [[namespace-shadows-mathlib]]
