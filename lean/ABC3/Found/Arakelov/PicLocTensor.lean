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

/-! ## ★出典の紐付け(`.src`) -/

def localizedTensorEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所化とテンソルの交換)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
