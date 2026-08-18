import ABC3.Found.Arakelov.PicTildeStalk
import Mathlib.RingTheory.Localization.BaseChange

/-!
# Arakelov (B1) 第 79 ブロック —— **局所化とテンソルの交換**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★比較射の茎の中身

    M_S ⊗_{R_S} N_S  ≅  (M ⊗_R N)_S

★これが `tilde` の比較射が茎で同型であることの**中身**である。

## ★★機構 —— mathlib の在庫 2 つの合成

| 段 | 在庫 |
|---|---|
| `M_S ⊗_{R_S} N_S ≅ M_S ⊗_R N_S` | ★`IsLocalization.moduleTensorEquiv` |
| `M_S ⊗_R N_S ≅ (M ⊗_R N)_S` | ★★**instance** `IsLocalizedModule S (TensorProduct.map f g)` |

★★★2 段目は mathlib が **instance として持っている**(2026-08-18 実測)
——`Localization/BaseChange.lean` 行 113。

★係数を `R` から `Localization S` へ上げるのは
`LinearEquiv.extendScalarsOfIsLocalization`。

## ★★本ブロックで取れるもの

| 定義 | 内容 |
|---|---|
| `localizedTensorEquiv` | ★★★★**`M_S ⊗_{R_S} N_S ≃ₗ (M ⊗_R N)_S`** |
| `localizedTensorEquiv_mk` | ★純テンソルでの値 |
-/

universe u

namespace ABC3.Found.Arakelov

open TensorProduct IsLocalization

variable (R : Type u) [CommRing R] (S : Submonoid R)
  (M N : Type u) [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N]

/-- ★★★★**局所化はテンソルと交換する** —— `M_S ⊗_{R_S} N_S ≅ (M ⊗_R N)_S`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `tilde` の比較射が茎で同型であることの中身である。 -/
noncomputable def localizedTensorEquiv :
    LocalizedModule S M ⊗[Localization S] LocalizedModule S N
      ≃ₗ[Localization S] LocalizedModule S (M ⊗[R] N) :=
  IsLocalization.moduleTensorEquiv S (Localization S) (LocalizedModule S M)
      (LocalizedModule S N)
    ≪≫ₗ (IsLocalizedModule.iso S (TensorProduct.map (LocalizedModule.mkLinearMap S M)
        (LocalizedModule.mkLinearMap S N))).symm.extendScalarsOfIsLocalization S (Localization S)

/-- ★★**分母 1 の純テンソルでの値**。 -/
@[simp] theorem localizedTensorEquiv_mk_one (m : M) (n : N) :
    localizedTensorEquiv R S M N
        (LocalizedModule.mk m 1 ⊗ₜ[Localization S] LocalizedModule.mk n 1)
      = LocalizedModule.mk (m ⊗ₜ[R] n) 1 := by
  have h : (TensorProduct.map (LocalizedModule.mkLinearMap S M)
      (LocalizedModule.mkLinearMap S N)) (m ⊗ₜ[R] n)
      = LocalizedModule.mk m 1 ⊗ₜ[R] LocalizedModule.mk n 1 := rfl
  simp only [localizedTensorEquiv, IsLocalization.moduleTensorEquiv,
    TensorProduct.equivOfCompatibleSMul, LinearEquiv.trans_apply,
    LinearEquiv.extendScalarsOfIsLocalization_apply, LinearEquiv.coe_mk,
    TensorProduct.mapOfCompatibleSMul_tmul]
  rw [show (LocalizedModule.mk m 1 ⊗ₜ[R] LocalizedModule.mk n 1) = _ from h.symm]
  exact IsLocalizedModule.iso_symm_apply _ _ _

/-- ★`a • mk m a = mk m 1`。 -/
theorem smul_mk_self (m : M) (a : S) :
    (a : R) • LocalizedModule.mk m a = LocalizedModule.mk m 1 := by
  rw [LocalizedModule.smul'_mk]
  exact LocalizedModule.mk_cancel a m

/-- ★**`S` の元による作用は単射**——局所化では `S` が可逆に働く。 -/
theorem smul_injective (a : S) :
    Function.Injective (fun x : LocalizedModule S M => (a : R) • x) := by
  have h := IsLocalizedModule.map_units (LocalizedModule.mkLinearMap S M) a
  rw [Module.End.isUnit_iff] at h
  exact h.1

/-- ★`R` の作用は左の因子を通る。 -/
theorem smul_tmul_left (x : LocalizedModule S M) (y : LocalizedModule S N) (r : R) :
    r • (x ⊗ₜ[Localization S] y) = (r • x) ⊗ₜ[Localization S] y := by
  rw [← IsScalarTower.algebraMap_smul (Localization S) r x,
    ← IsScalarTower.algebraMap_smul (Localization S) r (x ⊗ₜ[Localization S] y),
    TensorProduct.smul_tmul']

/-- ★`R` の作用は右の因子を通る。 -/
theorem smul_tmul_right (x : LocalizedModule S M) (y : LocalizedModule S N) (r : R) :
    r • (x ⊗ₜ[Localization S] y) = x ⊗ₜ[Localization S] (r • y) := by
  rw [← IsScalarTower.algebraMap_smul (Localization S) r y,
    ← IsScalarTower.algebraMap_smul (Localization S) r (x ⊗ₜ[Localization S] y),
    TensorProduct.tmul_smul]

/-- ★★★★**分母つきの値** —— `mk m a ⊗ mk n b ↦ mk (m ⊗ n) (a·b)`。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが前層射の構成で**局所分数条件が積で閉じる**ことの中身である。 -/
theorem localizedTensorEquiv_mk (m : M) (n : N) (a b : S) :
    localizedTensorEquiv R S M N
        (LocalizedModule.mk m a ⊗ₜ[Localization S] LocalizedModule.mk n b)
      = LocalizedModule.mk (m ⊗ₜ[R] n) (a * b) := by
  refine smul_injective R S (M ⊗[R] N) (a * b) ?_
  show ((a * b : S) : R) • _ = ((a * b : S) : R) • _
  rw [smul_mk_self]
  have hmul : ((a * b : S) : R) = (a : R) * (b : R) := rfl
  have hL : ((a * b : S) : R)
      • (LocalizedModule.mk m a ⊗ₜ[Localization S] LocalizedModule.mk n b)
      = LocalizedModule.mk m 1 ⊗ₜ[Localization S] LocalizedModule.mk n 1 := by
    rw [hmul, mul_smul, smul_tmul_right, smul_mk_self, smul_tmul_left, smul_mk_self]
  have hlin : ((a * b : S) : R) • localizedTensorEquiv R S M N
        (LocalizedModule.mk m a ⊗ₜ[Localization S] LocalizedModule.mk n b)
      = localizedTensorEquiv R S M N (((a * b : S) : R)
        • (LocalizedModule.mk m a ⊗ₜ[Localization S] LocalizedModule.mk n b)) := by
    rw [← IsScalarTower.algebraMap_smul (Localization S) ((a * b : S) : R)
      (LocalizedModule.mk m a ⊗ₜ[Localization S] LocalizedModule.mk n b), map_smul,
      IsScalarTower.algebraMap_smul]
  rw [hlin, hL, localizedTensorEquiv_mk_one]

/-! ## ★出典の紐付け(`.src`) -/

def localizedTensorEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所化とテンソルの交換)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
