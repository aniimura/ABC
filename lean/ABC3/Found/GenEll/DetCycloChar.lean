/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.EllModuliGalois
import ABC3.Found.GaloisRep.FullImageWitness
import Mathlib.NumberTheory.Cyclotomic.CyclotomicCharacter
import ABC3.Meta.Claim

/-!
# `det ρ` は mathlib の円分指標そのものである（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.22。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

## ★★★★★★★★これは何か

`Found/GaloisRep/FullImageWitness.lean` の `det_cyclotomic_full`（★無条件）は

    σ ζ = ζ ^ (det ρ(σ) mod l^n)   （`ζ^(l^n) = 1`）

を与える。★mathlib の `cyclotomicCharacter L l : (L ≃+* L) →* ℤ_[l]ˣ` は
**同じ性質で特徴づけられる**（`cyclotomicCharacter.spec`）。
★★したがって両者は一致する——それが本ファイルの内容である。

## ★★★★なぜ要るか

`imageSurjectiveJ_of_containsSL2`（`§9-1187`、第 761）が受けている仮説は

    ∀ u : ℤ_l^×, ∃ σ, det ρ(σ) = u

であった。本ファイルにより、これは**楕円曲線を含まない**言明

    円分指標 `cyclotomicCharacter L̄ l : Gal(L̄/L) → ℤ_l^×` は全射

に書き換わる。★`Skeleton/GenEll/GaloisLocal.lean` の葉 4 がそれである。

## ★機構

`ζ` を原始 `l^n` 乗根に取ると `ζ^a = ζ^b`（`a, b < l^n`）から `a = b`、
すなわち `toZModPow n` が一致する。★`PadicInt.ext_of_toZModPow` で `ℤ_[l]` の元として一致。
-/

namespace ABC3.Found.GenEll

open ABC3.Found.GaloisRep ABC3.Interface.GaloisRep WeierstrassCurve Matrix
open scoped Classical

/-- ★★★★★★★★★★★★★★
**`det ρ(σ)` は mathlib の円分指標そのものである**——★**無条件**。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

★これで `imageSurjective_of_containsSL2` 欄に残る仮説が
**楕円曲線を含まない言明**（円分指標の全射性）に書き換わる。 -/
theorem det_galRep_eq_cyclotomicCharacter (E : SSCurve) (l : ℕ) [Fact l.Prime]
    (e : E.tate l ≃+ (Fin 2 → ℤ_[l])) (σ : E.alg ≃ₐ[E.fld] E.alg) :
    ((Matrix.GeneralLinearGroup.det (galRep E.W l e σ) : ℤ_[l]ˣ) : ℤ_[l])
      = ((cyclotomicCharacter E.alg l σ.toRingEquiv : ℤ_[l]ˣ) : ℤ_[l]) := by
  haveI : WeierstrassCurve.IsElliptic ((E.W.baseChange E.alg).toAffine) :=
    isElliptic_baseChange_affine E.W E.isEll
  refine PadicInt.ext_of_toZModPow.mp (fun n => ?_)
  obtain ⟨ζ, hζ⟩ := HasEnoughRootsOfUnity.exists_primitiveRoot E.alg (l ^ n)
  have hz : ζ ^ (l ^ n) = 1 := hζ.pow_eq_one
  have h1 := det_cyclotomic_full E.W l e σ n ζ hz
  have h2 := cyclotomicCharacter.spec (L := E.alg) l σ.toRingEquiv ζ hz
  haveI : NeZero (l ^ n) := ⟨pow_ne_zero _ (Fact.out (p := l.Prime)).ne_zero⟩
  have heq : ζ ^ ((PadicInt.toZModPow n
        ((galRep E.W l e σ : Matrix (Fin 2) (Fin 2) ℤ_[l]).det)).val)
      = ζ ^ (((cyclotomicCharacter E.alg l σ.toRingEquiv).val.toZModPow n).val) := by
    rw [← h1, ← h2]
    rfl
  have hval := hζ.pow_inj (ZMod.val_lt _) (ZMod.val_lt _) heq
  exact ZMod.val_injective _ hval

/-- ★★★★★★★★★★★★★★★★
**円分指標が全射なら、`SL₂` を含むことから全射性が出る**——★**無条件**。

原文 (GenEll p.22):
> Corollary 4.3. (Full Galois Actions for Degenerating Elliptic Curves)

★★これが `imageSurjective_of_containsSL2` 欄の**最終形**である。
仮説 `hcyc` は**楕円曲線を含まない**——原文の括弧
『`ℚ(ζ_{l^∞})/ℚ` は `l` で完全分岐するので `L/ℚ` と線型無関連であり、
円分指標が全射になる』そのものである。 -/
theorem imageSurjectiveJ_of_cyclotomic (E : SSCurve) (l : ℕ) [Fact l.Prime]
    (hcyc : ∀ u : ℤ_[l]ˣ, ∃ σ : E.alg ≃ₐ[E.fld] E.alg,
      cyclotomicCharacter E.alg l σ.toRingEquiv = u)
    (h : ImageContainsSL2J E l) : ImageSurjectiveJ E l := by
  refine imageSurjectiveJ_of_containsSL2 E l (fun _ e u => ?_) h
  obtain ⟨σ, hσ⟩ := hcyc u
  refine ⟨σ, Units.ext ?_⟩
  rw [det_galRep_eq_cyclotomicCharacter E l e σ, hσ]

def imageSurjectiveJ_of_cyclotomic.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(imageSurjective_of_containsSL2 欄の最終形——円分指標の全射性だけ)",
    sectionId := "genell-cor-4-3" }

def imageSurjectiveJ_of_cyclotomic.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("☆残るのは「l が L で不分岐なら円分指標 cyclotomicCharacter L̄ l が全射」だけである。" ++
       "★★これは**楕円曲線を含まない**言明であり、mathlib の " ++
       "IsCyclotomicExtension.autEquivPow と「l が L で不分岐なら cyclotomic (l^n) は " ++
       "L 上既約」から出る") 8 ]

/-! ## ★出典の紐付け(`.src`) -/

def det_galRep_eq_cyclotomicCharacter.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 22,
    item := "Corollary 4.3(det ρ は mathlib の円分指標そのもの。★無条件)",
    sectionId := "genell-cor-4-3" }

end ABC3.Found.GenEll
