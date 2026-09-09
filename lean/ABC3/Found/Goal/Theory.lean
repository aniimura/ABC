/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Meta.Claim

/-!
# 登記 —— `Goal`

★このディレクトリは**論文でも既存理論でもない**。
16 本の論文の頂点を 1 つの最終目標(abc 予想)へ繋ぐ**接続点**だけを置く層である。

★置いてよいもの: 最終目標の主張、各論文の頂点の型の引用、頂点から目標への還元。
★置いてはならないもの: 数学の中身。中身は論文のディレクトリか理論の層に属する。

繋がり具合の測定は `node tools/goal-chain.mjs`(計測器は `Check/GoalChain.lean`)。
-/

namespace ABC3.Found.Goal

/-- **Goal の登記**。 -/
def theory : ABC3.Meta.Theory :=
  { what := "最終目標(abc 予想)の主張と、各論文の頂点をそこへ繋ぐ還元だけを置く接続層"
    consumers := ["GenEll", "IUTchIII", "pGC"]
    mathlibStatus :=
      "mathlib に abc 予想の statement は 0 件(2026-09-09 実測、`abcConjecture` で 0 hit)。" ++
      "`Nat.primeFactors` はあるので根基 `rad` はそれで定義した。" }

end ABC3.Found.Goal
