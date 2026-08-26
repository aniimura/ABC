import ABC3.Found.GaloisRep.NeronMinimal
import Mathlib.RingTheory.DedekindDomain.FiniteAdeleRing
import Mathlib.NumberTheory.NumberField.Basic

/-!
# Galois (G7) 第 314 ブロック —— **★★★★★★★悪い素点は有限個**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ].

## ★★★★★★★到達点

> 数体 `L` 上の曲線 `W`(`Δ ≠ 0`)に対し、
> **Néron 指数が `0` でない素点は有限個**(`finite_bad_primes`)

★★★これで `ω_E`(第 311 の分数イデアル)の**台が有限**になり、
分数イデアルとして意味を持つ。

## ★★★★★★素点ごとの環は `valuationSubringAtPrime`

`Localization.AtPrime p` を `L` の中で見るには `Algebra` の配管が要る。
★mathlib の **`HeightOneSpectrum.valuationSubringAtPrime L p`**(`L` の部分環)を使うと、
`Algebra`・`IsFractionRing`・`IsDomain` が**すべて自動で付く**。
★★離散付値環であることは `IsLocalization.AtPrime.isDiscreteValuationRing_of_dedekind_domain`
から出る(この部分環は `p` での局所化そのものだから)。

## ★★★★★悪い素点の中身

`p` が悪い(指数 `≠ 0`)なら、次のどれかが起きている:

* `a₁, a₂, a₃, a₄, a₆` のどれかが `p` で**極を持つ**
* `Δ` または `Δ⁻¹` が `p` で極を持つ(= `v(Δ) ≠ 0`)

★★★★どれも `HeightOneSpectrum.Support`(極の集合)で、**mathlib が有限性を持っている**
(`Support.finite`)。★7 つの有限集合の合併で押さえられる。

## ★★★なぜ `Δ` と `Δ⁻¹` の両方か

`v(Δ) = 0` を言うのに「`Δ` が整」だけでは足りない(`v(Δ) > 0` かもしれない)。
★`Δ⁻¹` も整なら `v(Δ) ≤ 0` も出て、合わせて `v(Δ) = 0 < 12`——第 313 が効く。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `mem_primeSubring_iff` | ★★付値 `≤ 1` と部分環の一致 |
| `isIntegral_of_mem` | ★★★係数が入っていれば整モデル |
| `minimalExp_eq_zero_of_mem` | ★★★★★良い素点では指数 `0` |
| `finite_bad_primes` | ★★★★★★★**悪い素点は有限個** |
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField WeierstrassCurve

variable {L : Type} [Field L] [NumberField L]

/-- ★★素点 `p` での付値環(`L` の部分環)。 -/
noncomputable abbrev primeSubring (p : HeightOneSpectrum (𝓞 L)) : ValuationSubring L :=
  p.valuationSubringAtPrime L

/-- ★★付値 `≤ 1` と部分環の一致。 -/
theorem mem_primeSubring_iff (p : HeightOneSpectrum (𝓞 L)) (x : L) :
    x ∈ primeSubring p ↔ (p.valuation L) x ≤ 1 := by
  show x ∈ p.valuationSubringAtPrime L ↔ _
  rw [IsDedekindDomain.HeightOneSpectrum.valuationSubringAtPrime_eq_valuationSubring]
  exact Valuation.mem_valuationSubring_iff _ _

instance instIsDiscreteValuationRingPrimeSubring (p : HeightOneSpectrum (𝓞 L)) :
    IsDiscreteValuationRing (primeSubring p) :=
  IsLocalization.AtPrime.isDiscreteValuationRing_of_dedekind_domain _ p.ne_bot _

/-- ★★★係数が部分環に入っていれば整モデル。 -/
theorem isIntegral_of_mem (S : ValuationSubring L) (W : WeierstrassCurve L)
    (h1 : W.a₁ ∈ S) (h2 : W.a₂ ∈ S) (h3 : W.a₃ ∈ S) (h4 : W.a₄ ∈ S) (h6 : W.a₆ ∈ S) :
    WeierstrassCurve.IsIntegral (↥S) W := by
  refine ⟨⟨⟨⟨W.a₁, h1⟩, ⟨W.a₂, h2⟩, ⟨W.a₃, h3⟩, ⟨W.a₄, h4⟩, ⟨W.a₆, h6⟩⟩, ?_⟩⟩
  refine WeierstrassCurve.ext ?_ ?_ ?_ ?_ ?_ <;> rfl

set_option maxHeartbeats 1600000 in
/-- ★★★★★係数と `Δ`・`Δ⁻¹` が入っていれば Néron 指数は `0`。 -/
theorem minimalExp_eq_zero_of_mem (W : WeierstrassCurve L) (hΔ : W.Δ ≠ 0)
    (p : HeightOneSpectrum (𝓞 L))
    (h1 : W.a₁ ∈ primeSubring p) (h2 : W.a₂ ∈ primeSubring p) (h3 : W.a₃ ∈ primeSubring p)
    (h4 : W.a₄ ∈ primeSubring p) (h6 : W.a₆ ∈ primeSubring p)
    (hd : W.Δ ∈ primeSubring p) (hdi : W.Δ⁻¹ ∈ primeSubring p) :
    minimalExp (primeSubring p) W = 0 := by
  haveI : WeierstrassCurve.IsIntegral (primeSubring p) W :=
    isIntegral_of_mem _ W h1 h2 h3 h4 h6
  have hnn : 0 ≤ vAdd (tateDvrVal (primeSubring p) L) (Units.mk0 W.Δ hΔ) :=
    tateDvrVal_nonneg_of_mem _ ⟨⟨W.Δ, hd⟩, rfl⟩
  have hnn' : 0 ≤ vAdd (tateDvrVal (primeSubring p) L) (Units.mk0 W.Δ hΔ)⁻¹ := by
    refine tateDvrVal_nonneg_of_mem _ ⟨⟨W.Δ⁻¹, hdi⟩, ?_⟩
    show W.Δ⁻¹ = _
    rw [Units.val_inv_eq_inv_val]
    rfl
  rw [vAdd_inv] at hnn'
  exact minimalExp_eq_zero_of_vAdd_Delta_lt W hΔ (by omega)

set_option maxHeartbeats 1600000 in
/-- ★★★★★★★**Néron 指数が `0` でない素点は有限個**。

原文 (GenEll p.15):
> Definition 3.3. We shall refer to the positive integer vK (qE ) ∈ Z&gt;0 as the local height of E [or EK ]. -/
theorem finite_bad_primes (W : WeierstrassCurve L) (hΔ : W.Δ ≠ 0) :
    {p : HeightOneSpectrum (𝓞 L) | minimalExp (primeSubring p) W ≠ 0}.Finite := by
  refine Set.Finite.subset
    ((((IsDedekindDomain.HeightOneSpectrum.Support.finite (𝓞 L) W.a₁).union
        (IsDedekindDomain.HeightOneSpectrum.Support.finite (𝓞 L) W.a₂)).union
      ((IsDedekindDomain.HeightOneSpectrum.Support.finite (𝓞 L) W.a₃).union
        (IsDedekindDomain.HeightOneSpectrum.Support.finite (𝓞 L) W.a₄))).union
      (((IsDedekindDomain.HeightOneSpectrum.Support.finite (𝓞 L) W.a₆).union
        (IsDedekindDomain.HeightOneSpectrum.Support.finite (𝓞 L) W.Δ)).union
        (IsDedekindDomain.HeightOneSpectrum.Support.finite (𝓞 L) W.Δ⁻¹))) ?_
  intro p hp
  by_contra hnot
  simp only [Set.mem_union, not_or] at hnot
  obtain ⟨⟨⟨n1, n2⟩, ⟨n3, n4⟩⟩, ⟨n6, nd⟩, ndi⟩ := hnot
  refine hp ?_
  refine minimalExp_eq_zero_of_mem W hΔ p ?_ ?_ ?_ ?_ ?_ ?_ ?_
  all_goals
    rw [mem_primeSubring_iff]
    first
      | exact not_lt.1 n1
      | exact not_lt.1 n2
      | exact not_lt.1 n3
      | exact not_lt.1 n4
      | exact not_lt.1 n6
      | exact not_lt.1 nd
      | exact not_lt.1 ndi

/-! ## ★出典の紐付け(`.src`) -/

def finite_bad_primes.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Definition 3.3(Néron 指数が 0 でない素点は有限個)",
    sectionId := "genell-def-3-3" }

end ABC3.Found.GaloisRep
