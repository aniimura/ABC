/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Meta.Claim

/-!
# 理論の登記 —— `SixExp`

★この木は**論文ではなく理論**である。`lean/ABC3/<bucket>/SixExp/` に置かれた宣言は、
原典の項目ではなく、それを支える**既存理論**を積んだものである。

## ★★依存は包含ではない

この理論を消費する論文の下へ移してはならない。
`GenEll` が mathlib を使っていても mathlib は `GenEll/` の下に無いのと同じで、
**使う側と使われる側は兄弟として並ぶ**。

## ☆台帳を Lean で持つ理由

`Meta/Claim.lean` の冒頭に書いてあるとおり、本プロジェクトは台帳を JSON ではなく
**Lean の宣言そのもの**として持つ——ノードとファイルが別々に生まれないので、
被覆の乖離が構造的に起こりえない。
★本ファイルはその規約に従って、理論の所属を Lean 側に置く
(2026-09-03 に一度 `ResearchPaper/theories.json` へ書いたが、規約違反だったので戻した)。
-/

namespace ABC3.Found.SixExp

/-- **SixExp の登記**。 -/
def theory : ABC3.Meta.Theory :=
  { what := "6 指数定理まわりの解析(Schneider–Lang 型の超越性)"
    consumers := ["FrdI"]
    mathlibStatus := "mathlib に超越性の道具はほぼ無い。" }

end ABC3.Found.SixExp
