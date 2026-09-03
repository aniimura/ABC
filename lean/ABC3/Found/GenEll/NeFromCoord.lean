/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.GalPointCoord
import ABC3.Found.GaloisRep.TateUnipotent
import ABC3.Meta.Claim

/-!
# 第 1301 ブロック —— **`h1` は「座標が基礎体に無い `l`-捩れ点」から出る**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★これは何か——`h1` の完成形

第 1298（座標が基礎体に無ければ動かされる）と
第 1270（`E[l]` で非自明なら `T_l E` でも `mod l` 非自明）を繋ぐ。

☆入力は「`l`-捩れ点で、`x` 座標が `K₀` に無いものがある」だけである。
★Tate 曲線ではそれが **`E[l] ⊄ E(K₀)`**、すなわち `l ∤ v(Q)` から出る（第 1279）。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Found.GaloisRep
open ABC3.Meta

open scoped Classical

variable {K₀ M : Type} [Field K₀] [Field M] [Algebra K₀ M]

/-- ★★★★★★★★★★★★★★★★
**座標が基礎体に無い `l`-捩れ点があれば `h1` が出る**——★**無条件**（第 1301）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1298 で `σ` を取り、第 1270 で `T_l E` に上げる。 -/
theorem exists_galTate_ne_of_coord_not_mem [IsGalois K₀ M] [IsAlgClosed M] [CharZero M]
    (W : WeierstrassCurve K₀) [WeierstrassCurve.IsElliptic (W.baseChange M).toAffine]
    (l : ℕ) [Fact l.Prime]
    (P : (W.baseChange M).toAffine.Point) (hP : l • P = 0)
    (h : (pointCoords P).1 ∉ Set.range (algebraMap K₀ M)) :
    ∃ σ : M ≃ₐ[K₀] M, ∃ x : tateModule (W.baseChange M) l,
      ∀ u : tateModule (W.baseChange M) l, galTate W l σ x ≠ x + l • u := by
  obtain ⟨σ, hσ⟩ := exists_galPoint_ne_of_coord_not_mem W P h
  exact ⟨σ, exists_galTate_ne_of_galPoint W l σ P hP hσ⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_galTate_ne_of_coord_not_mem.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(座標が基礎体に無い l-捩れ点があれば h1 が出る。★無条件)",
    sectionId := "genell-thm-3-8" }

def exists_galTate_ne_of_coord_not_mem.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_galPoint_ne_of_coord_not_mem(第 1298、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.exists_galPoint_ne_of_coord_not_mem") 1,
    .citation "[ABC3]" "exists_galTate_ne_of_galPoint(第 1270、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_galTate_ne_of_galPoint") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1301）**——`h1` の完成形である。" ++
       "☆Tate 曲線では入力が **`E[l] ⊄ E(K₀)`**、すなわち `l ∤ v(Q)` から出る（第 1279）。") 2 ]

end ABC3.Found.GenEll
