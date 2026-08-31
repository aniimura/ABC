/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.EllModuliGalois
import ABC3.Found.GenEll.GalRepContinuity
import Mathlib.FieldTheory.Galois.Profinite
import ABC3.Meta.Claim

/-!
# `galRep` の像は閉——連続性だけに帰着する（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これは何か

`Found/GenEll/Sl2Padic.lean` の `Lemma 3.1, (iv)` は `J` が**閉部分群**であることを要求する。
`Found/GenEll/EllModuliGalois.lean` の `imageContainsSL2J_of_alpha`（`§9-1188`、第 762）は
それを仮説 `hclosed` として受けていた。★本ファイルはそれを

    `galRep` が連続であること

**だけ**に帰着する。

## ★★★★mathlib の測定（2026-08-31、第 765）

☆着手前は「`Gal(L̄/L)` が profinite であることを自分で作る必要がある」と見ていたが、
**mathlib にすべて揃っている**（実測）:

* `IsGalois E.fld E.alg` —— instance で出る
* `CompactSpace (E.alg ≃ₐ[E.fld] E.alg)` —— `Mathlib/FieldTheory/Galois/Profinite.lean`
* `T2Space (GL (Fin 2) ℤ_[l])` —— instance で出る

★したがって「コンパクトの連続像は（Hausdorff で）閉」の一行で済む。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Interface.GaloisRep WeierstrassCurve Matrix
open scoped Classical

/-- ★★★★★★★★★★**`galRep` が連続なら像は閉部分群**——★**無条件**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★`Gal(L̄/L)` は Krull 位相で**コンパクト**（`Mathlib/FieldTheory/Galois/Profinite.lean`）、
`GL₂(ℤ_l)` は Hausdorff なので、コンパクトの連続像は閉である。 -/
theorem galRep_range_isClosed_of_continuous (E : SSCurve) (l : ℕ) [Fact l.Prime]
    (e : E.tate l ≃+ (Fin 2 → ℤ_[l]))
    (hcont : Continuous (galRep E.W l e)) :
    IsClosed (((galRep E.W l e).range : Subgroup (GL (Fin 2) ℤ_[l])) :
      Set (GL (Fin 2) ℤ_[l])) := by
  have hset : (((galRep E.W l e).range : Subgroup (GL (Fin 2) ℤ_[l])) :
      Set (GL (Fin 2) ℤ_[l])) = Set.range (galRep E.W l e) := by
    ext x
    simp [MonoidHom.mem_range]
  rw [hset, ← Set.image_univ]
  exact (isCompact_univ.image hcont).isClosed

/-- ★★★★★★★★★★★★
**`imageContainsSL2_of_torsionExt` 欄——`galRep` の連続性と `α` だけに帰着**。 -/
theorem imageContainsSL2J_of_alpha_of_continuous (E : SSCurve) (l : ℕ) [Fact l.Prime]
    (hl5 : 5 ≤ l)
    (hcont : ∀ e : E.tate l ≃+ (Fin 2 → ℤ_[l]), Continuous (galRep E.W l e))
    (halpha : ∀ e : E.tate l ≃+ (Fin 2 → ℤ_[l]),
      (Matrix.SpecialLinearGroup.toGL (upper (1 : ZMod l)) : GL (Fin 2) (ZMod l))
        ∈ ((galRep E.W l e).range).map (glRedPadic l))
    (hno : ¬ HasLCyclicJ E l) :
    ImageContainsSL2J E l :=
  imageContainsSL2J_of_alpha E l hl5
    (fun e => galRep_range_isClosed_of_continuous E l e (hcont e)) halpha hno

/-! ## ★★★★★★★★★★★★★★★★★★★★★★到達点——仮説は `α` だけ -/

/-- ★★★★★★★★★★★★★★**`galRep` の像は閉部分群**——★**無条件**。 -/
theorem galRep_range_isClosed (E : SSCurve) (l : ℕ) [Fact l.Prime]
    (e : E.tate l ≃+ (Fin 2 → ℤ_[l])) :
    IsClosed (((galRep E.W l e).range : Subgroup (GL (Fin 2) ℤ_[l])) :
      Set (GL (Fin 2) ℤ_[l])) :=
  galRep_range_isClosed_of_continuous E l e (galRep_continuous' E.W l e)

/-- ★★★★★★★★★★★★★★★★★★★★★★
**`imageContainsSL2_of_torsionExt` 欄——残る仮説は `α` だけ**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★★位相の側（像の閉性）は `§9-1198`（第 772）で閉じたので、
残るのは**局所理論の行列表示**——`α = (1 1 / 0 1)` が mod `l` 像に入ること——だけである。 -/
theorem imageContainsSL2J_of_alpha' (E : SSCurve) (l : ℕ) [Fact l.Prime] (hl5 : 5 ≤ l)
    (halpha : ∀ e : E.tate l ≃+ (Fin 2 → ℤ_[l]),
      (Matrix.SpecialLinearGroup.toGL (upper (1 : ZMod l)) : GL (Fin 2) (ZMod l))
        ∈ ((galRep E.W l e).range).map (glRedPadic l))
    (hno : ¬ HasLCyclicJ E l) :
    ImageContainsSL2J E l :=
  imageContainsSL2J_of_alpha E l hl5 (fun e => galRep_range_isClosed E l e) halpha hno

/-! ## ★出典の紐付け(`.src`) -/

def imageContainsSL2J_of_alpha'.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(imageContainsSL2 欄——残る仮説は α だけ)",
    sectionId := "genell-thm-3-8" }

def galRep_range_isClosed_of_continuous.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(galRep が連続なら像は閉部分群。★無条件)",
    sectionId := "genell-thm-3-8" }

def imageContainsSL2J_of_alpha_of_continuous.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(imageContainsSL2 欄——galRep の連続性と α だけに帰着)",
    sectionId := "genell-thm-3-8" }

def imageContainsSL2J_of_alpha_of_continuous.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆残るのは (i) galRep の連続性(E[l^n] への作用が有限拡大を経由すること)と " ++
       "(ii) 局所理論の行列表示(α が mod l 像に入る)の 2 本である") 12 ]

end ABC3.Found.GenEll
