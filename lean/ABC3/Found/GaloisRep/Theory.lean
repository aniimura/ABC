/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Meta.Claim

/-!
# 理論の登記 —— `GaloisRep`

★この木は**論文ではなく理論**である。`lean/ABC3/<bucket>/GaloisRep/` に置かれた宣言は、
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

namespace ABC3.Found.GaloisRep

/-- **GaloisRep の登記**。

★[GenEll] Lemma 3.5 の 483 ブロックのうち大半がここに積まれた。依存グラフ上は `Lemma 3.5` の 1 節点にしか見えていない(第 1450 の測定)。 -/
def theory : ABC3.Meta.Theory :=
  { what := "楕円曲線の Galois 表現・Tate 一意化・Vélu の同種・半安定性・最小判別式"
    consumers := ["GenEll"]
    mathlibStatus := "Tate 曲線・Vélu の公式・モジュラー多項式はいずれも mathlib に無い(2026-09-02 実測、ResearchPaper/mathlib-gap.json)。★形式群は Mathlib/RingTheory/FormalGroup/Basic.lean に定義 147 行だけある。" }

end ABC3.Found.GaloisRep
