import ABC3.Found.Arakelov.PicAwayLin

/-!
# Arakelov (B1) 第 126 ブロック —— **同型を `𝒪(D f)` 係数へ移す**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★局所自明性の最後の材料

第 100 ブロックは `M_f ≅ R_f` を **`Localization (powers f)` 係数**で与える。
★局所自明性は **`𝒪(D f)` 係数**で要る。

★★環同型 `awayRingEquiv`(第 124)で移すが、
`extendScalarsOfIsLocalization` は使えない(係数環が**同型だが別物**)ので
**手で `LinearEquiv` を組む**。

★★★`left_inv` / `right_inv` は `show` で β 簡約してから
`apply_symm_apply` / `symm_apply_apply` を当てる(実測)。

## ★★本ブロックで取れるもの

| 定義 | 内容 |
|---|---|
| `awayEquivScalar` | ★★★★★**`M_f ≃ₗ[𝒪(D f)] 𝒪(D f)`** |

## ★★★これで第 103 ブロック(生成元の乗法は全単射)が `𝒪(D f)` の側で当たる
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (M : ModuleCat.{u} (R : Type u)) (f : (R : Type u))

/-- ★★★★★**`M_f ≃ₗ[𝒪(D f)] 𝒪(D f)`**——同型を `𝒪(D f)` 係数へ移す。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★係数環が「同型だが別物」なので `extendScalarsOfIsLocalization` は使えない
——手で組む。 -/
noncomputable def awayEquivScalar
    (e : LocalizedModule (Submonoid.powers f) M
      ≃ₗ[Localization (Submonoid.powers f)] Localization (Submonoid.powers f)) :
    LocalizedModule (Submonoid.powers f) M
      ≃ₗ[(Γ(Spec R, PrimeSpectrum.basicOpen f) : Type u)]
        (Γ(Spec R, PrimeSpectrum.basicOpen f) : Type u) where
  toFun := fun x => (awayRingEquiv R f).symm (e x)
  invFun := fun c => e.symm (awayRingEquiv R f c)
  left_inv := fun x => by
    show e.symm ((awayRingEquiv R f) ((awayRingEquiv R f).symm (e x))) = x
    rw [AlgEquiv.apply_symm_apply, LinearEquiv.symm_apply_apply]
  right_inv := fun c => by
    show (awayRingEquiv R f).symm (e (e.symm ((awayRingEquiv R f) c))) = c
    rw [LinearEquiv.apply_symm_apply, AlgEquiv.symm_apply_apply]
  map_add' := fun x y => by simp
  map_smul' := fun c x => by
    show (awayRingEquiv R f).symm (e ((awayRingEquiv R f c) • x)) = _
    rw [map_smul]
    show (awayRingEquiv R f).symm ((awayRingEquiv R f c) * (e x)) = c * _
    rw [map_mul, AlgEquiv.symm_apply_apply]

/-! ## ★出典の紐付け(`.src`) -/

def awayEquivScalar.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——同型を 𝒪(D f) 係数へ移す)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
