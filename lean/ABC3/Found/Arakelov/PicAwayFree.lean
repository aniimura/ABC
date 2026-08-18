import ABC3.Found.Arakelov.PicTildeTensor
import Mathlib.RingTheory.Localization.Free

/-!
# Arakelov (B1) 第 76 ブロック —— **可逆加群は基本開集合の上で自由**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★素点から基本開集合へ**広げる**

`equivPicRing` の 8 段(§9-75)の第 1 段である。

★可逆加群 `M` は各素点 `p` で `M_p ≅ R_p`(mathlib の instance)。
★★しかし層の言葉に載せるには**基本開集合** `D(r)` の上で自由でなければならない。

★★★**mathlib に広げる補題がある**(2026-08-18 実測):

    Module.FinitePresentation.exists_free_localizedModule_powers
      (`RingTheory/Localization/Free.lean` 行 78)

「`M` が有限表示で `M_S` が `R_S` 上自由なら、ある `r ∈ S` で `M_r` は既に `R_r` 上自由」。

## ★★機構

| 段 | 在庫 |
|---|---|
| `M` 可逆 ⟹ 有限表示 | ★`Module.finitePresentation_of_projective` |
| `M_p` は可逆かつ自由 | ★mathlib instance |
| 素点 → 基本開集合 | ★★★`exists_free_localizedModule_powers` |
| 自由 ⟺ 環と線型同型 | ★`Module.Invertible.free_iff_linearEquiv` |

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `exists_away_linearEquiv` | ★★★★**各素点の近傍に `M_r ≅ R_r` となる `D(r)`** |
| `exists_away_linearEquiv_pair` | ★★2 つの可逆加群を**同時に**自由にする `r` |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite

variable (R : Type u) [CommRing R]

/-- ★★★★**可逆加群は各素点の近傍の基本開集合の上で自由である**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★これが `equivPicRing` の第 1 段である——素点での自由性を
**基本開集合へ広げる**。 -/
theorem exists_away_linearEquiv (M : Type u) [AddCommGroup M] [Module R M]
    [Module.Invertible R M] (p : PrimeSpectrum R) :
    ∃ r : R, r ∉ p.asIdeal ∧
      Nonempty (LocalizedModule.Away r M ≃ₗ[Localization (Submonoid.powers r)]
        Localization (Submonoid.powers r)) := by
  haveI : Module.FinitePresentation R M := Module.finitePresentation_of_projective R M
  obtain ⟨r, hr, hfree, -⟩ :=
    Module.FinitePresentation.exists_free_localizedModule_powers
      p.asIdeal.primeCompl (LocalizedModule.mkLinearMap p.asIdeal.primeCompl M)
      (Localization p.asIdeal.primeCompl)
  exact ⟨r, hr, Module.Invertible.free_iff_linearEquiv.mp hfree⟩

/-- ★★**2 つの可逆加群を同時に自由にする基本開集合が取れる**。

★`r · s` を取れば良い——`D(r·s) = D(r) ⊓ D(s)` であり、素イデアルの補元は積で閉じる。 -/
theorem exists_away_linearEquiv_pair (M N : Type u) [AddCommGroup M] [Module R M]
    [Module.Invertible R M] [AddCommGroup N] [Module R N] [Module.Invertible R N]
    (p : PrimeSpectrum R) :
    ∃ r s : R, r ∉ p.asIdeal ∧ s ∉ p.asIdeal ∧
      Nonempty (LocalizedModule.Away r M ≃ₗ[Localization (Submonoid.powers r)]
        Localization (Submonoid.powers r)) ∧
      Nonempty (LocalizedModule.Away s N ≃ₗ[Localization (Submonoid.powers s)]
        Localization (Submonoid.powers s)) := by
  obtain ⟨r, hr, hM⟩ := exists_away_linearEquiv R M p
  obtain ⟨s, hs, hN⟩ := exists_away_linearEquiv R N p
  exact ⟨r, s, hr, hs, hM, hN⟩

/-! ## ★出典の紐付け(`.src`) -/

def exists_away_linearEquiv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——可逆加群は基本開集合の上で自由)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
