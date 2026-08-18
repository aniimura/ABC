import ABC3.Found.Arakelov.PicStalkIso

/-!
# Arakelov (B1) 第 78 ブロック —— **`tilde M` の茎は `M_p`**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★比較射の左辺の茎

第 77 ブロックで「茎で同型なら同型」が建った。★比較射

    tildeTensorCmp : tilde (M ⊗_R N) ⟶ tensorModules (tilde M) (tilde N)

の**左辺の茎**は mathlib が既に持っている:

    instance : IsLocalizedModule x.asIdeal.primeCompl (tilde.toStalk M x).hom

★★すなわち `(tilde M)` の茎は `M_p` である(局所化加群として)。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `tildeStalkLocalized` | ★`tilde M` の茎は `M_p`(mathlib の在庫を名付ける) |
| `tildeStalkEquiv` | ★★**`M_p ≃ₗ[R] (tilde M) の茎`** |

## ★★★残り —— 右辺の茎

★右辺は `tensorModules`(= 前層テンソルの層化)である。
★★**層化は茎を変えない**(mathlib `stalkFunctor_map_unit_toSheafify_isIso`)ので、
前層テンソルの茎を求めれば良い。

★★★そこは基本開集合が**共終**であることを使う:

    茎 = colim_{f ∉ p} (M_f ⊗_{R_f} N_f) ≅ colim_{f ∉ p} (M ⊗_R N)_f = (M ⊗_R N)_p

★`M_f ⊗_{R_f} N_f ≅ (M ⊗_R N)_f` は mathlib(`moduleTensorEquiv` + `rTensor`)。
★★共終性は mathlib の `exists_mem_germ_eq_of_isBasis` /
`stalkFunctor_map_injective_of_isBasis` で扱える(2026-08-18 実測)。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (M : ModuleCat.{u} (R : Type u))

/-- ★**`tilde M` の茎は `M_p` である**——mathlib の在庫に名を付ける。 -/
theorem tildeStalkLocalized (x : PrimeSpectrum.Top (R : Type u)) :
    IsLocalizedModule x.asIdeal.primeCompl (tilde.toStalk M x).hom := inferInstance

/-- ★★**`M_p ≃ₗ[R] (tilde M) の茎`**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが比較射の**左辺の茎**である。 -/
noncomputable def tildeStalkEquiv (x : PrimeSpectrum.Top (R : Type u)) :
    LocalizedModule x.asIdeal.primeCompl M ≃ₗ[R] ((tilde M).presheaf.stalk x) :=
  IsLocalizedModule.iso x.asIdeal.primeCompl (tilde.toStalk M x).hom

/-! ## ★出典の紐付け(`.src`) -/

def tildeStalkEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——tilde の茎は局所化)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
