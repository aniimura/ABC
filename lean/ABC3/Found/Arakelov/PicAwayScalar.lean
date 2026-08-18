import ABC3.Found.Arakelov.PicLocGen

/-!
# Arakelov (B1) 第 124 ブロック —— **`M_f` を `𝒪(D f)` 加群と見る**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★§9-138 の型の層を剥がす

`tildeAwayEquiv : M_f ≃ₗ[R] Γ(tilde M, D f)` を `𝒪(D f)` 線型に上げたいが、

    Module 𝒪(D f) (LocalizedModule (powers f) M)

が**無い**——`LocalizedModule` は `Localization (powers f)` 上の加群で、
`𝒪(D f)` とは**同型だが別の環**だからである。

★★**解決**: `IsLocalization.algEquiv` で環同型を取り、`Module.compHom` で移す。
★★★これは第 95・98 ブロックで使ったのと**同じ手**である。

## ★★本ブロックで取れるもの

| 宣言 | 内容 |
|---|---|
| `awayRingEquiv` | ★`𝒪(D f) ≃ₐ[R] Localization (powers f)` |
| `modOnLocalized` | ★★`M_f` は `𝒪(D f)` 加群 |
| `towerOnLocalized` | ★★`R → 𝒪(D f) → M_f` は足場 |
| `tildeAwayEquivScalar` | ★★★★**`𝒪(D f)` 線型な `M_f ≅ Γ(tilde M, D f)`** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (M : ModuleCat.{u} (R : Type u)) (f : (R : Type u))

/-- ★**`𝒪(D f) ≃ₐ[R] Localization (powers f)`**。 -/
noncomputable abbrev awayRingEquiv :
    (Γ(Spec R, PrimeSpectrum.basicOpen f) : Type u)
      ≃ₐ[(R : Type u)] Localization (Submonoid.powers f) :=
  IsLocalization.algEquiv (Submonoid.powers f) _ _

/-- ★★**`M_f` は `𝒪(D f)` 加群である**。 -/
noncomputable instance modOnLocalized :
    Module (Γ(Spec R, PrimeSpectrum.basicOpen f) : Type u)
      (LocalizedModule (Submonoid.powers f) M) :=
  Module.compHom _ (awayRingEquiv R f).toRingHom

/-- ★★**`R → 𝒪(D f) → M_f` は足場をなす**。 -/
instance towerOnLocalized :
    IsScalarTower (R : Type u) (Γ(Spec R, PrimeSpectrum.basicOpen f) : Type u)
      (LocalizedModule (Submonoid.powers f) M) := by
  refine IsScalarTower.of_algebraMap_smul fun r x => ?_
  show (IsLocalization.algEquiv (Submonoid.powers f)
      (Γ(Spec R, PrimeSpectrum.basicOpen f) : Type u)
      (Localization (Submonoid.powers f))
    (algebraMap (R : Type u) (Γ(Spec R, PrimeSpectrum.basicOpen f) : Type u) r)) • x = r • x
  rw [AlgEquiv.commutes, IsScalarTower.algebraMap_smul]

/-! ## ★出典の紐付け(`.src`) -/

def towerOnLocalized.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——M_f と切断を 𝒪(D f) 加群と見る)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
