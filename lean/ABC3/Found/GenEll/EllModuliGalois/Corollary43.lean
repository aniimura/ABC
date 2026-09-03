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
import ABC3.Found.GenEll.EllModuliGalois.Theorem38

/-!
# EllModuliGalois —— `[GenEll] Corollary 4.3` の分

☆もとの 1 枚を**条なしの項目ごと**に割ったものである（第 1458、案 a）。
★「1 ファイル = 1 ノード」を回復するための分割で、中身は動かしていない。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Interface.GaloisRep WeierstrassCurve Matrix
open Matrix.SpecialLinearGroup
open scoped MatrixGroups Classical

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

def ImageSurjectiveJ.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(ImageSurjective 欄——galRep が全射)",
    sectionId := "genell-cor-4-3" }

end ABC3.Found.GenEll
