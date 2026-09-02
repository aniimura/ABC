/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.VeluLatticeElliptic
import ABC3.Found.GenEll.PointVariableChange
import ABC3.Meta.Claim

/-!
# 第 1332 ブロック —— **`Q` の倍数の座標集合＝代表系の座標集合**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★これは何か——第 1330 を「点の集合」の語彙に直す

第 1330 は代表系 `T`（＝`{k·z₀}`）で書かれている。
★原文（と `veluQuotientFull` の消費側）は
`{pointCoords (k • Q) : 1 ≤ k < l}` で書く。

☆`Q = Φ(z₀)` かつ `Φ` は加法準同型（`uniformHom`、在庫）なので
`k • Q = Φ(k·z₀)` であり、`k·z₀ ∉ Λ`（`1 ≤ k < l`）では
`Φ` の座標は `(℘(k z₀), ℘′(k z₀)/2)` である（`uniformMap_of_notMem`、在庫）。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve ABC3.Meta

open scoped Classical

/-- ★★★★★★**`Φ` の座標は格子点の座標**——★**無条件**（第 1332）。 -/
theorem pointCoords_uniformMap (P : PeriodPair) (hΔ : latticeDisc P ≠ 0) {z : ℂ}
    (hz : z ∉ P.lattice) :
    pointCoords (uniformMap P hΔ z) = (latticePointX P z, latticePointY P z) := by
  rw [uniformMap_of_notMem P hΔ hz]
  rfl

/-- ★★★★★★★★★★★★
**`Q` の倍数の座標集合は代表系の座標集合**——★**無条件**（第 1332）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

☆`k • Q = Φ(k z₀)` と `uniformMap_of_notMem` を `Finset.image_image` で束ねるだけである。 -/
theorem image_pointCoords_nsmul_eq (P : PeriodPair) (hΔ : latticeDisc P ≠ 0)
    (z₀ : ℂ) (l : ℕ)
    (hnot : ∀ k ∈ (Finset.range l).erase 0, ((k : ℂ)) * z₀ ∉ P.lattice) :
    ((Finset.range l).erase 0).image
        (fun k : ℕ => pointCoords (k • uniformMap P hΔ z₀))
      = ((Finset.range l).erase 0).image
        (fun k : ℕ => (latticePointX P ((k : ℂ) * z₀), latticePointY P ((k : ℂ) * z₀))) := by
  refine Finset.image_congr ?_
  intro k hk
  show pointCoords (k • uniformMap P hΔ z₀)
    = (latticePointX P ((k : ℂ) * z₀), latticePointY P ((k : ℂ) * z₀))
  have hk' : k ∈ (Finset.range l).erase 0 := Finset.mem_coe.1 hk
  have hsm : k • uniformMap P hΔ z₀ = uniformMap P hΔ ((k : ℂ) * z₀) := by
    have h1 : k • (uniformHom P hΔ) z₀ = (uniformHom P hΔ) (k • z₀) := (map_nsmul _ _ _).symm
    have h2 : (k : ℂ) * z₀ = k • z₀ := by
      rw [nsmul_eq_mul]
    rw [h2]
    simpa using h1
  rw [hsm, pointCoords_uniformMap P hΔ (hnot k hk')]

/-! ## ★出典の紐付け(`.src`) -/

def pointCoords_uniformMap.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Φ の座標は格子点の座標。★無条件)",
    sectionId := "genell-lemma-3-5" }

def image_pointCoords_nsmul_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Q の倍数の座標集合は代表系の座標集合。★無条件)",
    sectionId := "genell-lemma-3-5" }

def image_pointCoords_nsmul_eq.needs : List ProofObligation :=
  [ .citation "[ABC3]" "uniformHom(在庫、Φ は加法準同型)"
      (.inProject "ABC3" "ABC3.Found.GenEll.uniformHom") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1332）**——第 1330 を「点の集合」の語彙に直す段である。" ++
       "☆これで `veluQuotientFull` の消費側（原文の `{k • Q}` の形）と繋がる。") 2 ]

end ABC3.Found.GenEll
