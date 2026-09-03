/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.DegInfLocal
import ABC3.Found.GaloisRep.AdicCompleteIntegers

/-!
# 第 964 ブロック —— **★★★★★★★★★★★★★★★★★★★★完備化の付値の橋**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か——(D3) の (b1)

本プロジェクトの局所の定理群は仮説

    `hp : ∀ x : L, v_{𝔪_R}(x) = v_p(x)`

を受けている。★これは `HeightOneSpectrum.valuation Lv 𝔪_R` の言葉であり、
mathlib の `valuedAdicCompletion_eq_valuation` は `Valued.v` の言葉である。
☆その 2 つを同一視する橋が要った。

★道は 3 段である:

1. **同値**——両方の単位球はどちらも `R = adicCompletionIntegers`。
   `Valuation.isEquiv_of_val_le_one` で同値になる。
2. **同値かつ両方全射なら等しい**——mathlib に見当たらなかったので内製した。
   ☆一意化元 `π`（`v π = exp(-1)`）と `ρ`（`w ρ = exp(-1)`）を取り、
   `w π = exp(-1)` を挿んでから、任意の `x` を `π^{-n}` と比べる。
3. `valuedAdicCompletion_eq_valuation`（mathlib）で `L` の側に戻す。

★★**これで `Lemma 3.5` の局所の定理群を数体の素点の完備化に
実際に当てられる**。

| 定理 | 内容 |
|---|---|
| `le_exp_neg_one_of_lt_one` | ★`a < 1` なら `a ≤ exp(-1)`（離散だから） |
| `exp_neg_one_zpow` | ★`exp(-1)^k = exp(-k)` |
| `eq_of_isEquiv_of_surjective` | ★★★★★★★★★★★★**同値＋全射 ⇒ 等しい** |
| `valuation_maximalIdeal_adicCompletion_eq` | ★★★★★★★★★★★★★★★★**2 つの付値は一致** |
| `valuation_algebraMap_adicCompletion` | ★★★★★★★★★★★★★★★★★★★★**`hp` そのもの** |
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField IsDedekindDomain.HeightOneSpectrum
open scoped WithZero

/-! ## ★`ℤᵐ⁰` の小道具 -/

/-- ★**`a < 1` なら `a ≤ exp(-1)`**——値が離散（`ℤ` の像）だから。 -/
theorem le_exp_neg_one_of_lt_one (a : WithZero (Multiplicative ℤ)) (h : a < 1) :
    a ≤ WithZero.exp (-1 : ℤ) := by
  rcases eq_or_ne a 0 with rfl | ha
  · exact zero_le
  · obtain ⟨k, rfl⟩ : ∃ k : ℤ, a = WithZero.exp k :=
      ⟨WithZero.log a, (WithZero.exp_log ha).symm⟩
    rw [WithZero.exp_le_exp]
    have h0 : (WithZero.exp k : WithZero (Multiplicative ℤ)) < WithZero.exp (0 : ℤ) := by
      simpa using h
    rw [WithZero.exp_lt_exp] at h0
    omega

/-- ★**`exp(-1)^k = exp(-k)`**。 -/
theorem exp_neg_one_zpow (k : ℤ) :
    (WithZero.exp (-1 : ℤ)) ^ k = WithZero.exp (-k) := by
  simp only [WithZero.exp]
  rw [← WithZero.coe_zpow]
  congr 1
  rw [← ofAdd_zsmul]
  congr 1
  simp

/-! ## ★★★★★★★★★★★★同値＋全射 ⇒ 等しい -/

/-- ★★★★★★★★★★★★**同値で両方全射な `ℤᵐ⁰`-付値は等しい**。

☆mathlib には見当たらなかった（2026-09-01 の探索）ので内製する。

★道は 2 段:

1. `v π = exp(-1)` なる `π` を取ると、`w π = exp(-1)` も成り立つ。
   ☆`w ρ = exp(-1)` なる `ρ` を使う——`v ρ < 1` なので `v ρ ≤ v π`、
   同値で `w ρ ≤ w π`、一方 `w π < 1` なので `w π ≤ exp(-1)`。
2. 任意の `x ≠ 0` は `v x = v (π^{-n})` となる `n` を持ち、
   `IsEquiv.eq_iff` で `w x = w (π^{-n})`。 -/
theorem eq_of_isEquiv_of_surjective {K : Type} [Field K]
    {v w : Valuation K (WithZero (Multiplicative ℤ))} (h : v.IsEquiv w)
    (hv : Function.Surjective v) (hw : Function.Surjective w) : v = w := by
  obtain ⟨π, hπ⟩ := hv (WithZero.exp (-1 : ℤ))
  obtain ⟨ρ, hρ⟩ := hw (WithZero.exp (-1 : ℤ))
  have hexplt : (WithZero.exp (-1 : ℤ) : WithZero (Multiplicative ℤ)) < 1 := by
    have h0 : (WithZero.exp (-1 : ℤ) : WithZero (Multiplicative ℤ))
        < WithZero.exp (0 : ℤ) := by
      rw [WithZero.exp_lt_exp]; omega
    simpa using h0
  have hvπlt : v π < 1 := by rw [hπ]; exact hexplt
  have hwρlt : w ρ < 1 := by rw [hρ]; exact hexplt
  have hwπlt : w π < 1 := by
    have h1 : v π < v 1 := by simpa using hvπlt
    simpa using (h.lt_iff_lt (x := π) (y := 1)).1 h1
  have hvρlt : v ρ < 1 := by
    have h1 : w ρ < w 1 := by simpa using hwρlt
    simpa using (h.lt_iff_lt (x := ρ) (y := 1)).2 h1
  have hle1 : v ρ ≤ v π := by rw [hπ]; exact le_exp_neg_one_of_lt_one _ hvρlt
  have hle2 : w ρ ≤ w π := (h.le_iff_le).1 hle1
  have hwπ : w π = WithZero.exp (-1 : ℤ) :=
    le_antisymm (le_exp_neg_one_of_lt_one _ hwπlt) (by rw [← hρ]; exact hle2)
  have hπ0 : π ≠ 0 := by
    intro hc
    rw [hc, map_zero] at hπ
    exact absurd hπ.symm (by simp [WithZero.exp])
  ext x
  rcases eq_or_ne x 0 with rfl | hx
  · simp
  · obtain ⟨n, hn⟩ : ∃ n : ℤ, v x = WithZero.exp n :=
      ⟨WithZero.log (v x), (WithZero.exp_log ((v.ne_zero_iff).2 hx)).symm⟩
    have hvp : v (π ^ (-n)) = WithZero.exp n := by
      rw [map_zpow₀, hπ, exp_neg_one_zpow, neg_neg]
    have hwx : w x = w (π ^ (-n)) := (h.eq_iff).1 (by rw [hn, hvp])
    rw [hwx, map_zpow₀, hwπ, exp_neg_one_zpow, neg_neg, hn]

/-! ## ★★★★★★★★★★★★★★★★完備化の 2 つの付値は一致する -/

/-- ★★★★★★★★★★★★★★★★**`𝔪_R` の付値と `Valued.v` は一致する**。

☆両方の単位球はどちらも `adicCompletionIntegers` であり、
どちらも全射である。 -/
theorem valuation_maximalIdeal_adicCompletion_eq (L : Type) [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) :
    (IsDiscreteValuationRing.maximalIdeal (p.adicCompletionIntegers L)).valuation
        (p.adicCompletion L)
      = (Valued.v : Valuation (p.adicCompletion L) (WithZero (Multiplicative ℤ))) := by
  refine eq_of_isEquiv_of_surjective ?_ ?_ ?_
  · refine Valuation.isEquiv_of_val_le_one (fun x => ?_)
    constructor
    · intro hx
      obtain ⟨a, ha⟩ := IsDiscreteValuationRing.exists_lift_of_le_one hx
      rw [← ha]
      exact ((mem_adicCompletionIntegers (𝓞 L) L p).1 a.2)
    · intro hx
      exact HeightOneSpectrum.valuation_le_one _
        (⟨x, (mem_adicCompletionIntegers (𝓞 L) L p).2 hx⟩ : p.adicCompletionIntegers L)
  · exact HeightOneSpectrum.valuation_surjective _ _
  · exact valuedAdicCompletion_surjective L p

/-- ★★★★★★★★★★★★★★★★★★★★**局所の定理群の仮説 `hp` そのもの**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 964）**——これが (D3) の (b1) である。
☆`minDeltaExp_eq_mul_of_tateParamR`（第 892）や `minDeltaExp_eq_mul_of_veluMu`（第 904）の
`hp` はこれで埋まる。 -/
theorem valuation_algebraMap_adicCompletion (L : Type) [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (x : L) :
    (HeightOneSpectrum.valuation (p.adicCompletion L)
        (IsDiscreteValuationRing.maximalIdeal (p.adicCompletionIntegers L)))
        (algebraMap L (p.adicCompletion L) x)
      = (HeightOneSpectrum.valuation L p) x := by
  rw [valuation_maximalIdeal_adicCompletion_eq L p]
  have h := valuedAdicCompletion_eq_valuation' (R := 𝓞 L) (K := L) p x
  convert h using 2
  rfl

/-! ## ★出典の紐付け(`.src`) -/

def eq_of_isEquiv_of_surjective.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(同値で両方全射な ℤᵐ⁰-付値は等しい。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valuation_maximalIdeal_adicCompletion_eq.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(𝔪_R の付値と Valued.v は一致する。★無条件)",
    sectionId := "genell-lemma-3-5" }

def valuation_algebraMap_adicCompletion.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(局所の定理群の仮説 hp そのもの。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
