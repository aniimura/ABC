/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Sl2Padic
import ABC3.Meta.Claim

/-!
# `SL₂ ⊆ J` かつ `det J = R^×` なら `J = GL₂`（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.22。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

## ★★★★★★★★これは何か

`Interface/GenEll/EllModuli.lean` の

    imageSurjective_of_containsSL2 : ∀ E l, Prime l → PrimeToRamification E l →
      ImageContainsSL2 E l → ImageSurjective E l

欄の**群論の側**である。★原文 p.22 は

> if `l` is any prime number that is **unramified in `L`**, then the image of the Galois
> representation … contains `SL₂(ℤ_l)` **if and only if** the Galois representation …
> is **surjective**

と書き、理由を括弧で『`ℚ(ζ_{l^∞})/ℚ` は `l` で完全分岐するので `L/ℚ` と線型無関連であり、
円分指標が全射になる』と述べている。

## ★★★分解

| 段 | 内容 | 状態 |
|---|---|---|
| 群論 | `SL₂ ⊆ J` ＋ `det J = R^×` ⟹ `J = GL₂` | ★**本ファイル**（無条件） |
| 数論 | `l` が `L` で不分岐なら `det ∘ ρ`（円分指標）は全射 | ☆残る |

★群論の側は 3 行である——`g` に対して `det h = det g` なる `h ∈ J` を取れば
`g·h⁻¹` の行列式は `1`、すなわち `SL₂` の元だから `J` に入り、`g = (g·h⁻¹)·h ∈ J`。
-/

namespace ABC3.Found.GenEll

open Matrix Matrix.SpecialLinearGroup
open scoped MatrixGroups

/-- ★★★★★★★★**`SL₂ ⊆ J` かつ行列式が全射なら `J` は全体**——★**無条件**。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

★`imageSurjective_of_containsSL2` 欄の**群論の側**である。 -/
theorem mem_of_sl2_of_det {R : Type*} [CommRing R] (J : Subgroup (GL (Fin 2) R))
    (hsl : ∀ g : SL(2, R), (toGL g : GL (Fin 2) R) ∈ J)
    (hdet : ∀ u : Rˣ, ∃ h : GL (Fin 2) R, h ∈ J ∧ Matrix.GeneralLinearGroup.det h = u)
    (g : GL (Fin 2) R) : g ∈ J := by
  obtain ⟨h, hhJ, hdeth⟩ := hdet (Matrix.GeneralLinearGroup.det g)
  have hd : Matrix.GeneralLinearGroup.det (g * h⁻¹) = 1 := by
    rw [map_mul, map_inv, hdeth, mul_inv_cancel]
  have hmat : ((g * h⁻¹ : GL (Fin 2) R) : Matrix (Fin 2) (Fin 2) R).det = 1 := by
    have := congrArg Units.val hd
    simpa using this
  have hgl : (toGL (⟨((g * h⁻¹ : GL (Fin 2) R) : Matrix (Fin 2) (Fin 2) R), hmat⟩ :
      SL(2, R)) : GL (Fin 2) R) = g * h⁻¹ := by
    ext i j
    rfl
  have hmem : (g * h⁻¹) ∈ J := hgl ▸ hsl _
  have := J.mul_mem hmem hhJ
  simpa using this

/-- ★★★★★★**表現の水準での形**——像が `SL₂` を含み行列式が全射なら全射。 -/
theorem surjective_of_sl2_of_det {R G : Type*} [CommRing R] [Group G]
    (ρ : G →* GL (Fin 2) R)
    (hsl : ∀ g : SL(2, R), (toGL g : GL (Fin 2) R) ∈ ρ.range)
    (hdet : ∀ u : Rˣ, ∃ σ : G, Matrix.GeneralLinearGroup.det (ρ σ) = u) :
    Function.Surjective ρ := by
  intro g
  have hmem : g ∈ ρ.range := by
    refine mem_of_sl2_of_det ρ.range hsl (fun u => ?_) g
    obtain ⟨σ, hσ⟩ := hdet u
    exact ⟨ρ σ, ⟨σ, rfl⟩, hσ⟩
  exact hmem

/-! ## ★出典の紐付け(`.src`) -/

def mem_of_sl2_of_det.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(SL₂ ⊆ J かつ行列式が全射なら J = GL₂。★無条件)",
    sectionId := "genell-cor-4-3" }

def surjective_of_sl2_of_det.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(imageSurjective_of_containsSL2 欄の群論の側。★無条件)",
    sectionId := "genell-cor-4-3" }

def surjective_of_sl2_of_det.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆残るのは数論の側——l が L で不分岐なら det ∘ ρ(円分指標)が全射であること。" ++
       "原文の括弧『ℚ(ζ_{l^∞})/ℚ は l で完全分岐するので L/ℚ と線型無関連』") 6 ]

end ABC3.Found.GenEll
