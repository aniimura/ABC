/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Meta.Claim
import Mathlib.Data.Real.Basic
import Mathlib.Data.Set.Finite.Basic
import ABC3.Interface.GenEll.HeightTheory.Definition11

/-!
# HeightTheory —— `[GenEll] Proposition 1.7` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Interface.GenEll

open ABC3.Meta

/-- **`Proposition 1.7` の被覆の設定**を受ける `Interface`。

原文 (GenEll p.9):
> be a generically finite morphism of normal, Z-proper, Z-flat schemes of dimension two.

★★**なぜ条件 (a)–(d) を 1 つの不透明な `Prop` にするのか。**

原文の条件は「`D_ℚ`, `E_ℚ` が reduced」「`D_ℚ = φ_ℚ^{-1}(E_ℚ)_red`」
「`φ_ℚ` が有限**エタール**」「分岐指数が `e` を割る」——**すべて scheme 論の語彙**であり、
`HeightTheoryData` は点と高さしか持たないので展開できない。

★**ここで条件を落とすことは許されない**——落とせば主張が強くなり、G5(主張の弱化禁止)の
裏返しで**偽の skeleton** になる。ゆえに「展開できないが落とさない」唯一の書き方として、
条件を **1 つの `Prop` フィールドとして持ち、定理の仮定に置く**。
★この手当ては `prop_1_7.needs` に `.implicitStep` として明記してある。 -/
structure CoveringSetup where
  /-- 被覆の側 `Y`。 -/
  DY : HeightTheoryData
  /-- 底の側 `Z`。 -/
  DZ : HeightTheoryData
  /-- `φ : Y → Z` が `ℚ̄` 値点に誘導する写像。 -/
  toPoint : DY.Point → DZ.Point
  /-- `D ⊆ Y`。 -/
  divY : DY.Divisor
  /-- `E ⊆ Z`。 -/
  divZ : DZ.Divisor
  /-- 原文「Let `e` be a positive integer」。 -/
  e : ℕ
  /-- `e` は正。 -/
  e_pos : 0 < e
  /-- 原文の条件 (a)(b)(c)(d) と「生成的有限射・正規・ℤ-固有・ℤ-平坦・次元 2」。
  ★**不透明な 1 つの述語として持つ**(上の docstring を参照)。 -/
  hyp : Prop
  /-- 原文「`φ_ℚ` restricts to a finite étale morphism `(U_Y)_ℚ → (U_Z)_ℚ`」の
  点レベルの帰結——`U_Y(ℚ̄)` は `U_Z(ℚ̄)` へ写る。 -/
  maps_compl : ∀ x ∈ DY.compl divY, toPoint x ∈ DZ.compl divZ

/-- ★Track B は何を作らねばならないか。 -/
def CoveringSetup.waiting : WaitingFor :=
  { what := "次元 2 の正規・ℤ-固有・ℤ-平坦スキームの生成的有限射と、その分岐指数・エタール性・(−)_red。★`hyp` を不透明な Prop のままにしているのは『展開できないが落とさない』ためであり、展開が Track B の仕事である"
    trackB := "★★2026-08-27 に更新 —— **2026-08-16 の見立て『scheme の言葉に持ち上げる前に環の言葉で書ける』は当たった**。Found/GenEll/DifferentKummer.lean(差積イデアルと Kummer 理論)と BDSlack.lean で Proposition 1.7 の (i)-(iv) を構成へ載せ替えた(第 423 ブロック)。★残るのは **局所から大域への組み立て**であって、scheme 論の持ち上げではない" }

def CoveringSetup.src : Source :=
  { paper := "GenEll", pdfPage := 9, item := "Proposition 1.7",
    sectionId := "genell-prop-1-7" }

end ABC3.Interface.GenEll
