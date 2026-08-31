/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.HtFaltJField
import Mathlib.FieldTheory.IntermediateField.Adjoin.Basic
import ABC3.Meta.Claim

/-!
# 2 つの埋め込みの合成体（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★★★★これは何か

`HtFaltJField.lean`（`§9-1167`、第 740）の `htFaltOf_congr_j_of_common` は
「共通の体 `L₃` へ上げられるなら」という形だった。本ファイルは **`L₃` を実際に作る**。

    L₁ →+* ℂ,  L₂ →+* ℂ  ⟹  L₃ := (L₁ の像) ⊔ (L₂ の像) ⊆ ℂ

## ★★★機構

* `AlgHom.fieldRange` で像を `IntermediateField ℚ ℂ` として取る
* `AlgEquiv.ofInjectiveField` で `L ≃ₐ[ℚ] (像)`——有限次元性はこれで移す
* `IntermediateField.finiteDimensional_sup`（mathlib の instance）で合成体も有限次元
* `IntermediateField.inclusion` で `L₁ →ₐ[ℚ] L₃`

☆`Algebra (𝓞 L₁) L₃`・`IsScalarTower (𝓞 L₁) L₁ L₃`・`IsScalarTower (𝓞 L₁) (𝓞 L₃) L₃`・
`Module.Finite (𝓞 L₁) (𝓞 L₃)` は **mathlib がすべて自動で出す**（測定済み）。
本ファイルが作るのは `Algebra L₁ L₃` と `IsScalarTower ℚ L₁ L₃` だけでよい。

## ★★★★★★★★到達点

    htFaltOf_congr_j_of_emb :
      E₁/L₁ と E₂/L₂ が半安定、ℂ の中で j が一致 ⟹ ht^Falt(E₁) = ht^Falt(E₂)

★★★これが `EllModuliData` の `faltingsHeight : ℂ → ℝ` の**well-defined 性**である。
-/

namespace ABC3.Found.GaloisRep

open NumberField IsDedekindDomain WeierstrassCurve

/-! ## ★★★★★埋め込みの像 -/

section Emb

variable (L : Type) [Field L] [NumberField L] (e : L →+* ℂ)

/-- ★★数体からの環準同型は自動的に `ℚ`-代数準同型である。 -/
noncomputable def ratAlgHom : L →ₐ[ℚ] ℂ :=
  { e with
    commutes' := fun r => congrArg (fun f : ℚ →+* ℂ => f r)
      (Subsingleton.elim (e.comp (algebraMap ℚ L)) (algebraMap ℚ ℂ)) }

@[simp] theorem ratAlgHom_apply (x : L) : ratAlgHom L e x = e x := rfl

/-- ★★★**埋め込みの像**——`ℂ` の中の中間体として。 -/
noncomputable def embField : IntermediateField ℚ ℂ := (ratAlgHom L e).fieldRange

/-- ★★★★`L` は像と `ℚ`-代数として同型。 -/
noncomputable def embEquiv : L ≃ₐ[ℚ] ↥(embField L e) :=
  AlgEquiv.ofInjectiveField (ratAlgHom L e)

@[simp] theorem embEquiv_coe (x : L) : ((embEquiv L e x : ↥(embField L e)) : ℂ) = e x := rfl

instance finiteDimensional_embField : FiniteDimensional ℚ ↥(embField L e) :=
  LinearEquiv.finiteDimensional (embEquiv L e).toLinearEquiv

instance numberField_embField : NumberField ↥(embField L e) where
  to_charZero := inferInstance
  to_finiteDimensional := inferInstance

end Emb

/-! ## ★★★★★★★★合成体 -/

section Two

variable (L₁ L₂ : Type) [Field L₁] [NumberField L₁] [Field L₂] [NumberField L₂]
  (e₁ : L₁ →+* ℂ) (e₂ : L₂ →+* ℂ)

/-- ★★★★★★**2 つの埋め込みの像が生成する合成体**。 -/
noncomputable def compositum : IntermediateField ℚ ℂ := embField L₁ e₁ ⊔ embField L₂ e₂

instance finiteDimensional_compositum :
    FiniteDimensional ℚ ↥(compositum L₁ L₂ e₁ e₂) :=
  IntermediateField.finiteDimensional_sup _ _

instance numberField_compositum : NumberField ↥(compositum L₁ L₂ e₁ e₂) where
  to_charZero := inferInstance
  to_finiteDimensional := inferInstance

/-- ★★★`L₁` から合成体への `ℚ`-代数準同型。 -/
noncomputable def toCompositum₁ : L₁ →ₐ[ℚ] ↥(compositum L₁ L₂ e₁ e₂) :=
  (IntermediateField.inclusion (le_sup_left : embField L₁ e₁ ≤ compositum L₁ L₂ e₁ e₂)).comp
    (embEquiv L₁ e₁).toAlgHom

/-- ★★★`L₂` から合成体への `ℚ`-代数準同型。 -/
noncomputable def toCompositum₂ : L₂ →ₐ[ℚ] ↥(compositum L₁ L₂ e₁ e₂) :=
  (IntermediateField.inclusion (le_sup_right : embField L₂ e₂ ≤ compositum L₁ L₂ e₁ e₂)).comp
    (embEquiv L₂ e₂).toAlgHom

@[simp] theorem toCompositum₁_coe (x : L₁) :
    ((toCompositum₁ L₁ L₂ e₁ e₂ x : ↥(compositum L₁ L₂ e₁ e₂)) : ℂ) = e₁ x := rfl

@[simp] theorem toCompositum₂_coe (x : L₂) :
    ((toCompositum₂ L₁ L₂ e₁ e₂ x : ↥(compositum L₁ L₂ e₁ e₂)) : ℂ) = e₂ x := rfl

noncomputable instance algebraCompositum₁ : Algebra L₁ ↥(compositum L₁ L₂ e₁ e₂) :=
  (toCompositum₁ L₁ L₂ e₁ e₂).toRingHom.toAlgebra

noncomputable instance algebraCompositum₂ : Algebra L₂ ↥(compositum L₁ L₂ e₁ e₂) :=
  (toCompositum₂ L₁ L₂ e₁ e₂).toRingHom.toAlgebra

instance isScalarTower_compositum₁ : IsScalarTower ℚ L₁ ↥(compositum L₁ L₂ e₁ e₂) :=
  IsScalarTower.of_algebraMap_eq (fun r => ((toCompositum₁ L₁ L₂ e₁ e₂).commutes r).symm)

instance isScalarTower_compositum₂ : IsScalarTower ℚ L₂ ↥(compositum L₁ L₂ e₁ e₂) :=
  IsScalarTower.of_algebraMap_eq (fun r => ((toCompositum₂ L₁ L₂ e₁ e₂).commutes r).symm)

theorem algebraMap_compositum₁ (x : L₁) :
    ((algebraMap L₁ ↥(compositum L₁ L₂ e₁ e₂) x : ↥(compositum L₁ L₂ e₁ e₂)) : ℂ) = e₁ x := rfl

theorem algebraMap_compositum₂ (x : L₂) :
    ((algebraMap L₂ ↥(compositum L₁ L₂ e₁ e₂) x : ↥(compositum L₁ L₂ e₁ e₂)) : ℂ) = e₂ x := rfl

/-! ## ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★到達点 -/

/-- ★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★★
**半安定な楕円曲線の `ht^Falt` は `ℂ` の中の `j` だけで決まる**——★**無条件**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

    E₁/L₁ 半安定、E₂/L₂ 半安定、e₁(j(E₁)) = e₂(j(E₂)) ⟹ ht^Falt(E₁) = ht^Falt(E₂)

★★★★これが `EllModuliData` の `faltingsHeight : EllClass → ℝ`（`EllClass := ℂ`、
`cls E = j(E)`）の **well-defined 性そのもの**である
（`ResearchPaper/ellmoduli-witness-status.json` の `designChoice`）。 -/
theorem htFaltOf_congr_j_of_emb (E₁ : WeierstrassCurve L₁) [E₁.IsElliptic]
    (E₂ : WeierstrassCurve L₂) [E₂.IsElliptic]
    (h₁ : ∀ p : HeightOneSpectrum (𝓞 L₁), SemistableAt p E₁)
    (h₂ : ∀ p : HeightOneSpectrum (𝓞 L₂), SemistableAt p E₂)
    (hj : e₁ E₁.j = e₂ E₂.j) :
    htFaltOf L₁ E₁ = htFaltOf L₂ E₂ := by
  refine htFaltOf_congr_j_of_common L₁ L₂ ↥(compositum L₁ L₂ e₁ e₂) E₁ E₂ h₁ h₂ ?_
  refine Subtype.ext ?_
  rw [algebraMap_compositum₁, algebraMap_compositum₂, hj]

end Two

/-! ## ★出典の紐付け(`.src`) -/

def compositum.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(2 つの埋め込みの像が生成する合成体)",
    sectionId := "genell-prop-3-4" }

def htFaltOf_congr_j_of_emb.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(半安定な曲線の ht^Falt は ℂ の中の j だけで決まる。★無条件)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GaloisRep
