/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Meta.Claim

/-!
# 理論の登記 —— `NumberField`

★この木は**論文ではなく理論**である。`lean/ABC3/<bucket>/NumberField/` に置かれた宣言は、
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

namespace ABC3.Found.NumberField

/-- **NumberField の登記**。

☆馴分岐(tame ramification)は `ResearchPaper/mathlib-gap.json` で 6/6 節点まで分解済み。 -/
def theory : ABC3.Meta.Theory :=
  { what := "数体の道具——判別式・different・分岐・素点の輸送"
    consumers := ["FrdI", "GenEll"]
    mathlibStatus := "`NumberField` の基礎は mathlib にある。ここは足りない補題の置き場である。" }

end ABC3.Found.NumberField
