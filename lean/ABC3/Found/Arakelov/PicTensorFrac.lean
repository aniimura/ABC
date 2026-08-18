import ABC3.Found.Arakelov.PicLocTensor

/-!
# Arakelov (B1) 第 81 ブロック —— **局所分数条件は積で閉じる**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.3。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

## ★★★★★前層射の構成の核

`tilde M` の切断は**明示的**である(2026-08-18 実測、`rfl` で確認):

    (tilde M)(U) = { f : Π x : U, M_x // isLocallyFraction f }

★★2 つの切断 `s`・`t` から点ごとにテンソルを取った

    (s ⊗ t)(x) := localizedTensorEquiv (s x ⊗ₜ t x)  ∈ (M ⊗_R N)_x

が**再び局所分数条件を満たす**——これが本ブロックである。

## ★★機構 —— mathlib の `add_mem'` の写経

    s x = mk r ⟨a⟩,  t x = mk r' ⟨b⟩   ⟹   (s ⊗ t)(x) = mk (r ⊗ₜ r') ⟨a·b⟩

★`sectionsSubmodule` の `add_mem'`(`sb • ra + sa • rb` / `sa * sb`)と**同じ形**である。
★★値の計算は第 80 ブロックの `localizedTensorEquiv_mk`。

## ★★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `tensorFun` | ★点ごとのテンソル |
| `tensorFun_frac` | ★★★★**局所分数条件は積で閉じる** |
-/

universe u

namespace ABC3.Found.Arakelov

open AlgebraicGeometry CategoryTheory Opposite TopologicalSpace TensorProduct StructureSheaf

variable (R : Type u) [CommRing R] (M N : Type u)
  [AddCommGroup M] [Module R M] [AddCommGroup N] [Module R N]
  (U : Opens (PrimeSpectrum.Top R))

/-- ★**点ごとのテンソル**。 -/
noncomputable def tensorFun
    (s : Π x : U, Localizations (R := R) M x.1)
    (t : Π x : U, Localizations (R := R) N x.1) :
    Π x : U, Localizations (R := R) (M ⊗[R] N) x.1 :=
  fun x => localizedTensorEquiv R x.1.asIdeal.primeCompl M N
    (s x ⊗ₜ[Localization x.1.asIdeal.primeCompl] t x)

/-- ★★★★**局所分数条件は積で閉じる**。

原文 (GenEll p.3):
> (i) We shall refer to as an arithmetic line bundle L = (L, | − |L) on X any

★★★`s x = mk r ⟨a⟩`・`t x = mk r' ⟨b⟩` なら `mk (r ⊗ₜ r') ⟨a·b⟩` である
——mathlib の `sectionsSubmodule.add_mem'` と同じ形。 -/
theorem tensorFun_frac
    {s : Π x : U, Localizations (R := R) M x.1}
    {t : Π x : U, Localizations (R := R) N x.1}
    (hs : (isLocallyFraction R M).pred s) (ht : (isLocallyFraction R N).pred t) :
    (isLocallyFraction R (M ⊗[R] N)).pred (tensorFun R M N U s t) := by
  intro x
  obtain ⟨Va, ma, ia, ra, sa, wa⟩ := hs x
  obtain ⟨Vb, mb, ib, rb, sb, wb⟩ := ht x
  refine ⟨Va ⊓ Vb, ⟨ma, mb⟩, Opens.infLELeft _ _ ≫ ia, ra ⊗ₜ[R] rb, sa * sb, fun y ↦ ?_⟩
  obtain ⟨hsay, hsa⟩ := wa ⟨y.1, y.2.1⟩
  obtain ⟨hsby, hsb⟩ := wb ⟨y.1, y.2.2⟩
  refine ⟨y.1.asIdeal.primeCompl.mul_mem hsay hsby, ?_⟩
  exact (congrArg (localizedTensorEquiv R _ M N)
      congr($hsa ⊗ₜ[Localization (y.1 : PrimeSpectrum R).asIdeal.primeCompl] $hsb)).trans
    (localizedTensorEquiv_mk R _ M N ra rb ⟨sa, hsay⟩ ⟨sb, hsby⟩)

/-! ## ★出典の紐付け(`.src`) -/

def tensorFun_frac.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 3,
    item := "Definition 1.1, (i)(層 B——局所分数条件は積で閉じる)",
    sectionId := "genell-def-1-1-i" }

end ABC3.Found.Arakelov
