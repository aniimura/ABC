/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.EllModuliObjects
import ABC3.Found.GaloisRep.GalRep
import ABC3.Found.GenEll.Sl2Padic
import ABC3.Found.GenEll.GLSurjective
import ABC3.Found.GenEll.Thm38Bridge
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

/-- ★★★★★**`ImageSurjective` 欄**——表現が全射であること。 -/
def ImageSurjectiveJ (E : SSCurve) (l : ℕ) : Prop :=
  ∀ (_ : Fact l.Prime) (e : E.tate l ≃+ (Fin 2 → ℤ_[l])),
    Function.Surjective (galRep E.W l e)

/-- ★★全射なら `SL₂` を含む（自明な向き）。 -/
theorem imageContainsSL2J_of_surjective (E : SSCurve) (l : ℕ)
    (h : ImageSurjectiveJ E l) : ImageContainsSL2J E l :=
  fun hl e g => h hl e (toGL g)

/-! ## ★★★★★★★★★★`imageSurjective_of_containsSL2` 欄の帰着 -/

/-- ★★★★★★★★★★★★
**行列式（円分指標）が全射なら、`SL₂` を含むことから全射性が出る**——★**無条件**。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

★★これが `imageSurjective_of_containsSL2` 欄の**全体**である。
群論の側は `Found/GenEll/GLSurjective.lean` の `surjective_of_sl2_of_det`（`§9-1186`）、
残るのは仮説 `hdet`——**`l` が `L` で不分岐なら円分指標が全射**——だけである。

☆`det ρ(σ)` が円分指標であること自体は
`Found/GaloisRep/FullImageWitness.lean` の `det_cyclotomic_full`（★無条件）で済んでいる。 -/
theorem imageSurjectiveJ_of_containsSL2 (E : SSCurve) (l : ℕ)
    (hdet : ∀ (_ : Fact l.Prime) (e : E.tate l ≃+ (Fin 2 → ℤ_[l])) (u : ℤ_[l]ˣ),
      ∃ σ : E.alg ≃ₐ[E.fld] E.alg,
        Matrix.GeneralLinearGroup.det (galRep E.W l e σ) = u)
    (h : ImageContainsSL2J E l) : ImageSurjectiveJ E l := by
  intro hl e
  exact surjective_of_sl2_of_det _ (fun g => h hl e g) (fun u => hdet hl e u)

def imageSurjectiveJ_of_containsSL2.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(imageSurjective_of_containsSL2 欄——円分指標の全射性に帰着)",
    sectionId := "genell-cor-4-3" }

def imageSurjectiveJ_of_containsSL2.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆残るのは「l が L で不分岐なら円分指標 Gal(L̄/L) → ℤ_l^× が全射」だけである。" ++
       "原文の括弧『ℚ(ζ_{l^∞})/ℚ は l で完全分岐するので L/ℚ と線型無関連』。" ++
       "★det ρ(σ) が円分指標であること自体は det_cyclotomic_full で済んでいる") 8 ]

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
    (halpha : ∀ e : E.tate l ≃+ (Fin 2 → ℤ_[l]),
      (toGL (upper (1 : ZMod l)) : GL (Fin 2) (ZMod l))
        ∈ ((galRep E.W l e).range).map (glRedPadic l))
    (hno : ¬ HasLCyclicJ E l) :
    ImageContainsSL2J E l := by
  intro _ e g
  have hnoe : ¬ ∃ e : E.tate l ≃+ (Fin 2 → ℤ_[l]),
      ∃ v : Fin 2 → ZMod l, v ≠ 0 ∧
        ∀ M ∈ (↑(((galRep E.W l e).range).map (glRedPadic l)) :
            Set (GL (Fin 2) (ZMod l))),
          ∃ c : ZMod l, (M : Matrix (Fin 2) (Fin 2) (ZMod l)).mulVec v = c • v :=
    fun h => hno (fun _ => h)
  have hno' : ¬ ∃ v : Fin 2 → ZMod l, v ≠ 0 ∧
      ∀ M ∈ (↑(((galRep E.W l e).range).map (glRedPadic l)) :
          Set (GL (Fin 2) (ZMod l))),
        ∃ c : ZMod l, (M : Matrix (Fin 2) (Fin 2) (ZMod l)).mulVec v = c • v :=
    fun h => hnoe ⟨e, h⟩
  exact sl2_of_alpha_of_no_stable_line l hl5 _ (hclosed e) (halpha e) hno' g

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

def ImageSurjectiveJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(ImageSurjective 欄——galRep が全射)",
    sectionId := "genell-cor-4-3" }

end ABC3.Found.GenEll
