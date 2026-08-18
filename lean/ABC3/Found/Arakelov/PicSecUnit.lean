import ABC3.Found.Arakelov.PicToOpenRes

/-!
# Arakelov (B1) 第 119 ブロック —— **`powers g` は `D(h·g)` の切断に可逆に作用する**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★制限の両立を `IsLocalizedModule.ext` で言うための前提

`tildeAwayEquiv`(第 85)どうしが制限で対応することを言うには、
**`powers g` が `Γ(tilde M, D(h·g))` に可逆に作用する**ことが要る。

★★mathlib が**そのまま持っている**(2026-08-24 実測):

    Scheme.Modules.isUnit_algebraMap_end_of_le_basicOpen (f) (hf : U ≤ basicOpen f) :
      IsUnit (algebraMap R (End R Γ(M, U)) f)

★★★`D(h·g) = D(h) ⊓ D(g) ≤ D(g)` なのでそのまま当たる。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `isUnit_end_section` | ★`g` は `Γ(tilde M, D(h·g))` に可逆に作用 |
| `isUnit_end_powers_section` | ★★★**`powers g` が可逆に作用** |

## ★★★次

これで `IsLocalizedModule.ext` が使え、
`restriction ∘ tildeAwayEquiv_g = tildeAwayEquiv_{h·g} ∘ (局所化)` が出る。
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace

variable (R : CommRingCat.{u}) (M : ModuleCat.{u} (R : Type u)) (g h : (R : Type u))

/-- ★**`g` は `Γ(tilde M, D(h·g))` に可逆に作用する**。 -/
theorem isUnit_end_section : IsUnit (algebraMap (R : Type u)
    (Module.End (R : Type u)
      (Γ(tilde M, PrimeSpectrum.basicOpen (h * g)) : Type u)) g) := by
  refine Scheme.Modules.isUnit_algebraMap_end_of_le_basicOpen g ?_
  rw [PrimeSpectrum.basicOpen_mul]
  exact inf_le_right

/-- ★★★**`powers g` は `Γ(tilde M, D(h·g))` に可逆に作用する**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★これで `IsLocalizedModule.ext` が使える。 -/
theorem isUnit_end_powers_section (s : Submonoid.powers g) :
    IsUnit (algebraMap (R : Type u)
      (Module.End (R : Type u)
        (Γ(tilde M, PrimeSpectrum.basicOpen (h * g)) : Type u)) (s : (R : Type u))) := by
  obtain ⟨n, hn⟩ := s.2
  rw [← hn, map_pow]
  exact (isUnit_end_section R M g h).pow n

/-! ## ★出典の紐付け(`.src`) -/

def isUnit_end_powers_section.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——powers g は切断に可逆に作用する)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
