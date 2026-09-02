/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluLatticeElliptic
import ABC3.Meta.Claim

/-!
# 第 1331 ブロック —— **Vélu の商の楕円性は体の埋め込みで降りる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★これは何か——`ℂ` から数体へ降ろす一行

第 1330 で **`ℂ` の上では Vélu の商が楕円**であることが取れた。
★楕円性は `Δ ≠ 0`（体の上）であり、`Δ` は体の準同型で移る
（`veluQuotientFull_map`、在庫）ので、**埋め込みの先で `Δ ≠ 0` なら元でも `Δ ≠ 0`** である。

☆これで「数体上の Vélu の商は楕円」が `ℂ` の計算に帰着する。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve ABC3.Meta

open scoped Classical

/-- ★★★★★★★★★★★★
**Vélu の商の楕円性は体の埋め込みで降りる**——★**無条件**（第 1331）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`Δ` は体の準同型で移り、体の準同型は単射だからである。 -/
theorem isElliptic_veluQuotientFull_of_map {F A : Type*} [Field F] [Field A] (f : F →+* A)
    (W : WeierstrassCurve F) (S : Finset (F × F))
    (h : (veluQuotientFull (W.map f) (S.image (fun Q => (f Q.1, f Q.2)))).IsElliptic) :
    (veluQuotientFull W S).IsElliptic := by
  have hΔ : ((veluQuotientFull W S).map f).Δ ≠ 0 := by
    rw [veluQuotientFull_map]
    exact h.isUnit.ne_zero
  rw [WeierstrassCurve.map_Δ] at hΔ
  have h0 : (veluQuotientFull W S).Δ ≠ 0 := by
    intro hc
    rw [hc, map_zero] at hΔ
    exact hΔ rfl
  exact ⟨isUnit_iff_ne_zero.2 h0⟩

/-! ## ★出典の紐付け(`.src`) -/

def isElliptic_veluQuotientFull_of_map.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の商の楕円性は体の埋め込みで降りる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isElliptic_veluQuotientFull_of_map.needs : List ProofObligation :=
  [ .citation "[ABC3]" "veluQuotientFull_map(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.veluQuotientFull_map") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1331）**——第 1330（`ℂ` の上の Vélu の定理）を" ++
       "数体へ降ろす一行である。☆これで「数体上の Vélu の商は楕円」が" ++
       "`ℂ` の計算に帰着する。") 2 ]

end ABC3.Found.GenEll
