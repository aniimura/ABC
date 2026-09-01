/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.BasisRealize
import ABC3.Found.GenEll.EllModuliGalois
import ABC3.Found.GenEll.GalRepClosed
import ABC3.Meta.Claim

/-!
# 第 1237 ブロック —— **`α` が `mod l` 像に入る（`halpha` の形）**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か

第 1236 は「ある基底で行列が `α`」を**行列の等式**として与える。
☆`imageContainsSL2J_of_alpha` が要る `halpha` は
`toGL (upper 1) ∈ (galRep の像).map (glRedPadic l)` という**群の言葉**なので、
`Units.ext` で繋ぐ。

★★★これで `AlphaBridge` の II 側は
**`σ` の `mod l` の冪単性と非自明性だけ**に帰着した。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Interface.GaloisRep WeierstrassCurve ABC3.Meta
open Matrix.SpecialLinearGroup
open scoped MatrixGroups Classical

/-- ★★★★★★★★★★★★★★★★★★★★
**`α` が `mod l` 像に入る（`halpha` の形）**——★**無条件**（第 1237）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆仮説は `σ` の `T_l E` 上の作用が `mod l` で**冪単かつ非自明**であることだけ。

★★★これが `imageContainsSL2J_of_alpha`（在庫）の `halpha` そのものである。 -/
theorem alpha_mem_map_of_galTate (E : SSCurve) (l : ℕ) [Fact l.Prime]
    (e : E.tate l ≃+ (Fin 2 → ℤ_[l])) (σ : E.alg ≃ₐ[E.fld] E.alg)
    (h2 : ∀ x : E.tate l, ∃ u : E.tate l,
      galTate E.W l σ (galTate E.W l σ x) + x
        = galTate E.W l σ x + galTate E.W l σ x + l • u)
    (h1 : ∃ x : E.tate l, ∀ u : E.tate l, galTate E.W l σ x ≠ x + l • u) :
    ∃ e₀ : E.tate l ≃+ (Fin 2 → ℤ_[l]),
      (toGL (upper (1 : ZMod l)) : GL (Fin 2) (ZMod l))
        ∈ ((galRep E.W l e₀).range).map (glRedPadic l) := by
  obtain ⟨e₀, he₀⟩ := exists_basis_glRed_eq_alpha l E.W e σ h2 h1
  refine ⟨e₀, Subgroup.mem_map.2 ⟨galRep E.W l e₀ σ,
    MonoidHom.mem_range.2 ⟨σ, rfl⟩, Units.ext ?_⟩⟩
  rw [he₀]
  rfl

/-! ## ★★★★★★★★★★★★★★★★★★★★`ImageContainsSL2J` へ -/

/-- ★★★★★★★★★★★★★★★★★★★★
**像が `SL₂(ℤ_l)` を含む（`galTate` の言葉で）**——★（第 1246）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`imageContainsSL2J_of_alpha`（在庫）の `halpha` を第 1237 で埋めた形。

★★★これが `EllModuliWitness` の `imageContainsSL2J_torsionExt` が要る形である
——残る仮説は `hclosed`（`galRep` の像が閉部分群）と、
`σ` の `mod l` の冪単性・非自明性だけ。 -/
theorem imageContainsSL2J_of_galTate (E : SSCurve) (l : ℕ) [Fact l.Prime] (hl5 : 5 ≤ l)
    (hclosed : ∀ e : E.tate l ≃+ (Fin 2 → ℤ_[l]),
      IsClosed (((galRep E.W l e).range : Subgroup (GL (Fin 2) ℤ_[l])) :
        Set (GL (Fin 2) ℤ_[l])))
    (e : E.tate l ≃+ (Fin 2 → ℤ_[l])) (σ : E.alg ≃ₐ[E.fld] E.alg)
    (h2 : ∀ x : E.tate l, ∃ u : E.tate l,
      galTate E.W l σ (galTate E.W l σ x) + x
        = galTate E.W l σ x + galTate E.W l σ x + l • u)
    (h1 : ∃ x : E.tate l, ∀ u : E.tate l, galTate E.W l σ x ≠ x + l • u)
    (hno : ¬ HasLCyclicJ E l) :
    ImageContainsSL2J E l :=
  imageContainsSL2J_of_alpha E l hl5 hclosed (alpha_mem_map_of_galTate E l e σ h2 h1) hno

/-- ★★★★★★★★★★★★★★★★★★★★★★★★
**像が `SL₂(ℤ_l)` を含む（`hclosed` も落とした形）**——★（第 1249）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆`imageContainsSL2J_of_alpha'`（在庫）は像の閉性を**すでに落としている**
（`galRep` の連続性は `galRep_continuous'`、第 772・葉 5）。

★★★したがって残る仮説は
**`σ` の `mod l` の冪単性・非自明性と `¬ HasLCyclicJ` だけ**である
——`EllModuliWitness` の `imageContainsSL2J_torsionExt` に必要なものはこれで尽きる。 -/
theorem imageContainsSL2J_of_galTate' (E : SSCurve) (l : ℕ) [Fact l.Prime] (hl5 : 5 ≤ l)
    (e : E.tate l ≃+ (Fin 2 → ℤ_[l])) (σ : E.alg ≃ₐ[E.fld] E.alg)
    (h2 : ∀ x : E.tate l, ∃ u : E.tate l,
      galTate E.W l σ (galTate E.W l σ x) + x
        = galTate E.W l σ x + galTate E.W l σ x + l • u)
    (h1 : ∃ x : E.tate l, ∀ u : E.tate l, galTate E.W l σ x ≠ x + l • u)
    (hno : ¬ HasLCyclicJ E l) :
    ImageContainsSL2J E l :=
  imageContainsSL2J_of_alpha' E l hl5 (alpha_mem_map_of_galTate E l e σ h2 h1) hno

/-! ## ★出典の紐付け(`.src`) -/

def imageContainsSL2J_of_galTate'.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(像が SL₂(Z_l) を含む——hclosed も落とした形)",
    sectionId := "genell-thm-3-8" }

def imageContainsSL2J_of_galTate.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(像が SL₂(Z_l) を含む——galTate の言葉で。hclosed は仮説)",
    sectionId := "genell-thm-3-8" }

def alpha_mem_map_of_galTate.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(α が mod l 像に入る——halpha の形。★無条件)",
    sectionId := "genell-thm-3-8" }

def alpha_mem_map_of_galTate.needs : List ProofObligation :=
  [ .citation "[ABC3]" "exists_basis_glRed_eq_alpha(第 1236、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GaloisRep.exists_basis_glRed_eq_alpha") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1237）**——`imageContainsSL2J_of_alpha`（在庫）の" ++
       "`halpha` そのものである。☆仮説は `σ` の `T_l E` 上の作用が `mod l` で" ++
       "**冪単かつ非自明**であることだけ。" ++
       "★★★これで `AlphaBridge` の II 側はその 2 条件に帰着した。") 3 ]

end ABC3.Found.GenEll
