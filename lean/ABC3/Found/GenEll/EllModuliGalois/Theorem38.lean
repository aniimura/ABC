/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.GenEll.EllModuliObjects
import ABC3.Found.GaloisRep.GalRepBasis
import ABC3.Found.GaloisRep.GalRep
import ABC3.Found.GenEll.Sl2Padic
import ABC3.Found.GenEll.GLSurjective
import ABC3.Found.GenEll.Thm38Bridge
import ABC3.Meta.Claim

/-!
# EllModuliGalois —— `[GenEll] Theorem 3.8` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
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
  ∀ (_ : Fact l.Prime), ∃ e : E.tate l ≃+ (Fin 2 → ℤ_[l]),
    ∃ v : Fin 2 → ZMod l, v ≠ 0 ∧
      ∀ M ∈ (↑(((galRep E.W l e).range).map (glRedPadic l)) :
          Set (GL (Fin 2) (ZMod l))),
        ∃ c : ZMod l, (M : Matrix (Fin 2) (Fin 2) (ZMod l)).mulVec v = c • v

/-! ## ★★★★★★`ImageContainsSL2` 欄・`ImageSurjective` 欄 -/

/-- ★★★★★★**`ImageContainsSL2` 欄**——像が `SL₂(ℤ_l)` を含むこと。 -/
def ImageContainsSL2J (E : SSCurve) (l : ℕ) : Prop :=
  ∀ (_ : Fact l.Prime) (e : E.tate l ≃+ (Fin 2 → ℤ_[l])) (g : SL(2, ℤ_[l])),
    (toGL g : GL (Fin 2) ℤ_[l]) ∈ (galRep E.W l e).range

/-! ## ★★★★★★★★★★★★`imageContainsSL2_of_torsionExt` 欄の帰着 -/

/-- ★★★★★★★★★★★★★★
**`α` が mod `l` 像に入り、安定直線が無ければ `SL₂(ℤ_l)` を含む**——★**無条件**。

原文 (GenEll p.19):
> Then the image of the Galois representation Gal(Q[bb][bar]/L) → GL_2(Z[bb]_l) associated to

★★これが `imageContainsSL2_of_torsionExt` 欄の**群論の側の全体**である。
中身は `Found/GenEll/Thm38Bridge.lean` の `sl2_of_alpha_of_no_stable_line`
（＝ `Lemma 3.1, (iv)` ＋ `exists_nonUpper_of_no_stable_line`）。

☆残るのは仮説 `halpha`——**乗法還元の素点で局所高さが `l` で割れなければ
mod `l` 像が `α = (1 1 / 0 1)` を含む**——だけである。
★原文の『by the local theory (cf. the discussion preceding Lemma 3.2)』の中身であり、
`Found/GaloisRep/Lemma32Tate.lean`（Tate 一意化と `Lemma 3.2, (i)`）が素材である。

☆仮説 `hclosed` は像が閉部分群であること（profinite 群の連続像）。 -/
theorem imageContainsSL2J_of_alpha (E : SSCurve) (l : ℕ) [hlp : Fact l.Prime] (hl5 : 5 ≤ l)
    (hclosed : ∀ e : E.tate l ≃+ (Fin 2 → ℤ_[l]),
      IsClosed (((galRep E.W l e).range : Subgroup (GL (Fin 2) ℤ_[l])) :
        Set (GL (Fin 2) ℤ_[l])))
    (halpha : ∃ e₀ : E.tate l ≃+ (Fin 2 → ℤ_[l]),
      (toGL (upper (1 : ZMod l)) : GL (Fin 2) (ZMod l))
        ∈ ((galRep E.W l e₀).range).map (glRedPadic l))
    (hno : ¬ HasLCyclicJ E l) :
    ImageContainsSL2J E l := by
  intro _ e g
  obtain ⟨e₀, halpha₀⟩ := halpha
  have hnoe : ¬ ∃ e : E.tate l ≃+ (Fin 2 → ℤ_[l]),
      ∃ v : Fin 2 → ZMod l, v ≠ 0 ∧
        ∀ M ∈ (↑(((galRep E.W l e).range).map (glRedPadic l)) :
            Set (GL (Fin 2) (ZMod l))),
          ∃ c : ZMod l, (M : Matrix (Fin 2) (Fin 2) (ZMod l)).mulVec v = c • v :=
    fun h => hno (fun _ => h)
  have hno₀ : ¬ ∃ v : Fin 2 → ZMod l, v ≠ 0 ∧
      ∀ M ∈ (↑(((galRep E.W l e₀).range).map (glRedPadic l)) :
          Set (GL (Fin 2) (ZMod l))),
        ∃ c : ZMod l, (M : Matrix (Fin 2) (Fin 2) (ZMod l)).mulVec v = c • v :=
    fun h => hnoe ⟨e₀, h⟩
  -- ★基底 `e₀` で `SL₂ ⊆ 像` を出し、共役ですべての基底へ移す
  have hsl0 : ∀ g : SL(2, ℤ_[l]),
      (Matrix.SpecialLinearGroup.toGL g : GL (Fin 2) ℤ_[l]) ∈ (galRep E.W l e₀).range :=
    fun g => sl2_of_alpha_of_no_stable_line l hl5 _ (hclosed e₀) halpha₀ hno₀ g
  exact sl2_range_basis_transfer E.W l e₀ hsl0 e g

def imageContainsSL2J_of_alpha.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 19,
    item := "Theorem 3.8(imageContainsSL2_of_torsionExt 欄の群論の側——α と安定直線から SL₂)",
    sectionId := "genell-thm-3-8" }

def imageContainsSL2J_of_alpha.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆残るのは局所理論の行列表示——乗法還元の素点で局所高さが l で割れなければ " ++
       "mod l 像が α = (1 1 / 0 1) を含むこと。" ++
       "原文の『by the local theory (cf. the discussion preceding Lemma 3.2)』であり、" ++
       "Found/GaloisRep/Lemma32Tate.lean(Tate 一意化と Lemma 3.2, (i))が素材である") 10,
    .implicitStep "☆像が閉部分群であること(profinite 群の連続像)" 3 ]

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

end ABC3.Found.GenEll
