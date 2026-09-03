/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.GalTateUnipFull
import ABC3.Found.GenEll.GalPointRational
import ABC3.Meta.Claim

/-!
# 第 1300 ブロック —— **基礎体に `ζ_l` と位数 `l` の点があれば、どの `σ` も幂単**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★★★★★★★★★★★★★これは何か——`h2` の**完成形**

第 1296 は「`σ` が `ζ_l` と位数 `l` の点を固定する」ことを要求していた。
★**どちらも基礎体の側の条件に直せる**:

* `ζ_l ∈ K₀` なら `σ ζ = ζ`（`AlgEquiv.commutes`）
* 位数 `l` の点が `K₀` 有理なら `σ` はそれを固定する（第 1299）

☆したがって **`h2` は「基礎体に `ζ_l` と位数 `l` の点がある」だけから、
すべての `σ` について出る**。

★★★Tate 曲線ではこの 2 つがちょうど `μ_l ⊂ E(K₀)`（第 1297）である。
-/

namespace ABC3.Found.GenEll

open WeierstrassCurve WeierstrassCurve.Affine ABC3.Interface.GaloisRep ABC3.Found.GaloisRep
open ABC3.Meta

open scoped Classical

variable {K₀ M : Type} [Field K₀] [CharZero K₀] [Field M] [Algebra K₀ M] [IsAlgClosed M]

set_option maxHeartbeats 800000 in
/-- ★★★★★★★★★★★★★★★★★★★★
**基礎体に `ζ_l` と位数 `l` の点があれば、どの `σ` も `mod l` 幂単**——★**無条件**（第 1300）。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

☆第 1296 に第 1299（基礎体の点は固定される）と `AlgEquiv.commutes` を渡すだけである。 -/
theorem galTate_unipotent_of_rational (W : WeierstrassCurve K₀) [W.IsElliptic]
    [WeierstrassCurve.IsElliptic (W.baseChange M).toAffine]
    (l : ℕ) [Fact l.Prime]
    (e : tateModule (W.baseChange M) l ≃+ (Fin 2 → ℤ_[l])) (σ : M ≃ₐ[K₀] M)
    (ζ₀ : K₀) (hζ : IsPrimitiveRoot (algebraMap K₀ M ζ₀) l)
    (P₀ : W.toAffine.Point) (hP₀ : addOrderOf P₀ = l) :
    ∀ x : tateModule (W.baseChange M) l, ∃ u : tateModule (W.baseChange M) l,
      galTate W l σ (galTate W l σ x) + x
        = galTate W l σ x + galTate W l σ x + l • u := by
  haveI : (W.map (algebraMap K₀ M)).IsElliptic :=
    inferInstanceAs ((W.baseChange M).IsElliptic)
  refine galTate_unipotent_of_fixed_root_point W l e σ (algebraMap K₀ M ζ₀) hζ
    (σ.commutes ζ₀) (rhPoint (algebraMap K₀ M) W P₀) ?_ (galPoint_rhPoint_eq W σ P₀)
  exact (addOrderOf_rhPoint (algebraMap K₀ M) W P₀).trans hP₀

/-! ## ★出典の紐付け(`.src`) -/

def galTate_unipotent_of_rational.src : Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(基礎体に ζ_l と位数 l の点があればどの σ も mod l 幂単。★無条件)",
    sectionId := "genell-thm-3-8" }

def galTate_unipotent_of_rational.needs : List ProofObligation :=
  [ .citation "[ABC3]" "galTate_unipotent_of_fixed_root_point(第 1296、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.galTate_unipotent_of_fixed_root_point") 1,
    .citation "[ABC3]" "galPoint_rhPoint_eq(第 1299、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.galPoint_rhPoint_eq") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1300）**——`h2` の**完成形**である。" ++
       "☆Tate 曲線ではこの 2 つの条件がちょうど `μ_l ⊂ E(K₀)`（第 1297）である。") 3 ]

end ABC3.Found.GenEll
