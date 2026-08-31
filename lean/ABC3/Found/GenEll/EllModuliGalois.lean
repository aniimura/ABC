/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.EllModuliObjects
import ABC3.Found.GaloisRep.GalRep
import ABC3.Found.GenEll.Sl2Padic
import ABC3.Meta.Claim

/-!
# `EllModuliData` の Galois 側の語彙（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.19。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

## ★★★★★★★★これは何か

`EllModuliData` の残る 3 つの述語欄

| 欄 | 内容 |
|---|---|
| `HasLCyclic` | `E` が `l`-巡回部分群スキームをもつ |
| `ImageContainsSL2` | `Gal → GL₂(ℤ_l)` の像が `SL₂(ℤ_l)` を含む |
| `ImageSurjective` | その表現が全射 |

を、**本プロジェクトの `galRep`（`Found/GaloisRep/GalRep.lean`）で書く**。

## ★★★`HasLCyclic` の読み替え（`Interface/GenEll/EllModuli.lean` の測定）

★標数 0 では有限群スキームはすべてエタール（Cartier）であり、
エタール有限群スキームは有限 Galois 加群と圏同値である。したがって

> 「`l`-巡回部分群スキーム `H ⊆ E`」 ⟺ 「`E[l]` の中の `Gal`-安定な直線」

である。★★`Found/GenEll/Thm38Bridge.lean` の `exists_nonUpper_of_no_stable_line` は
まさにこの形（`∃ v ≠ 0, ∀ M ∈ 像, ∃ c, M·v = c·v`）を受け取る。

## ★量化の向き

★どの欄も「**すべての**両立する基底 `e` について」と取る。
`T_l E ≃ ℤ_l²` の同一視を変えると `ρ` は共役で変わるが、
「`SL₂` を含む」も「安定直線をもつ」も**共役で不変**なので、
`∀ e` と `∃ e` は同値である（その同値性の証明は本ファイルでは取らない）。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Interface.GaloisRep WeierstrassCurve Matrix
open Matrix.SpecialLinearGroup
open scoped MatrixGroups Classical

/-- `E` の定義体の代数閉包。 -/
abbrev SSCurve.alg (E : SSCurve) : Type := AlgebraicClosure E.fld

/-- `E` の `l` 進 Tate 加群（代数閉包の上で）。 -/
abbrev SSCurve.tate (E : SSCurve) (l : ℕ) : Type :=
  tateModule (E.W.baseChange E.alg) l

/-! ## ★★★★★★`HasLCyclic` 欄 -/

/-- ★★★★★★**`HasLCyclic` 欄**——`E[l]` の中に `Gal`-安定な直線があること。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★標数 0 では `l`-巡回部分群スキームと `Gal`-安定な直線は同じことである
（`Interface/GenEll/EllModuli.lean` の `HasLCyclic` の docstring の測定）。 -/
def HasLCyclicJ (E : SSCurve) (l : ℕ) : Prop :=
  ∀ (_ : Fact l.Prime) (e : E.tate l ≃+ (Fin 2 → ℤ_[l])),
    ∃ v : Fin 2 → ZMod l, v ≠ 0 ∧
      ∀ σ : E.alg ≃ₐ[E.fld] E.alg, ∃ c : ZMod l,
        ((glRedPadic l (galRep E.W l e σ) : GL (Fin 2) (ZMod l)) :
          Matrix (Fin 2) (Fin 2) (ZMod l)).mulVec v = c • v

/-! ## ★★★★★★`ImageContainsSL2` 欄・`ImageSurjective` 欄 -/

/-- ★★★★★★**`ImageContainsSL2` 欄**——像が `SL₂(ℤ_l)` を含むこと。 -/
def ImageContainsSL2J (E : SSCurve) (l : ℕ) : Prop :=
  ∀ (_ : Fact l.Prime) (e : E.tate l ≃+ (Fin 2 → ℤ_[l])) (g : SL(2, ℤ_[l])),
    (toGL g : GL (Fin 2) ℤ_[l]) ∈ (galRep E.W l e).range

/-- ★★★★★**`ImageSurjective` 欄**——表現が全射であること。 -/
def ImageSurjectiveJ (E : SSCurve) (l : ℕ) : Prop :=
  ∀ (_ : Fact l.Prime) (e : E.tate l ≃+ (Fin 2 → ℤ_[l])),
    Function.Surjective (galRep E.W l e)

/-- ★★全射なら `SL₂` を含む（自明な向き）。 -/
theorem imageContainsSL2J_of_surjective (E : SSCurve) (l : ℕ)
    (h : ImageSurjectiveJ E l) : ImageContainsSL2J E l :=
  fun hl e g => h hl e (toGL g)

/-! ## ★出典の紐付け(`.src`) -/

def HasLCyclicJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(HasLCyclic 欄——E[l] の中の Gal-安定な直線)",
    sectionId := "genell-thm-3-8" }

def HasLCyclicJ.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★標数 0 では有限群スキームはすべてエタール(Cartier)であり、" ++
       "エタール有限群スキームは有限 Galois 加群と圏同値なので、" ++
       "l-巡回部分群スキームと Gal-安定な直線は同じことである。" ++
       "☆その圏同値そのものは形式化していない——読み替えとして記録する") 4 ]

def ImageContainsSL2J.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(ImageContainsSL2 欄——galRep の像が SL₂(ℤ_l) を含む)",
    sectionId := "genell-thm-3-8" }

def ImageSurjectiveJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(ImageSurjective 欄——galRep が全射)",
    sectionId := "genell-cor-4-3" }

end ABC3.Found.GenEll
