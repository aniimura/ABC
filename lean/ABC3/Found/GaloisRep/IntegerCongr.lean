/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import Mathlib.RingTheory.DedekindDomain.FiniteAdeleRing
import Mathlib.NumberTheory.NumberField.Completion.FinitePlace
import ABC3.Meta.Claim

/-!
# 第 988 ブロック —— **★★★★★★★★★★★★★★★★`𝓞 L` は `p` を法として代表を尽くす**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★これは何か

第 982（分裂性の二者択一）は剰余体が**有限**であることを要求する。
☆だが `Finite (IsLocalRing.ResidueField (p.adicCompletionIntegers L))` は
mathlib に無い（第 983・985 で実測）。

★その構成の**核**が本ブロックである:

    `v_p(y) ≤ 1` なる `y ∈ L` に対し、`a ∈ 𝓞 L` を `v_p(y − a) < 1` に取れる

☆これがあれば、`L` の稠密性と合わせて `𝓞 L → ResidueField R` の全射性が出て、
`𝓞 L / p` の有限性（`Ideal.finiteQuotientOfFreeOfNeBot`）に帰着する。

★道は 4 段である:

| 段 | 中身 | 出どころ |
|---|---|---|
| 1 | `y` の分母に現れる素点は有限個 | `HeightOneSpectrum.Support.finite` |
| 2 | それらは `p` と異なる（`v_p(y) ≤ 1` だから） | 直接 |
| 3 | 中国剰余で `s ≡ 1 (mod p)`、`s ∈ q^{n_q}` を取る | `exists_forall_sub_mem_ideal` |
| 4 | `s·y` は整（全付値が `≤ 1`）で `y ≡ s·y (mod p)` | `mem_integers_of_valuation_le_one` |
-/

namespace ABC3.Found.GaloisRep

open IsDedekindDomain NumberField
open scoped Classical

/-- ★★★★★★★★★★★★★★★★**`v_p(y) ≤ 1` なら `y` は `𝓞 L` の元と `p` を法として合同**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 988）**——完備化の剰余体の有限性（mathlib の穴）を建てる核である。

☆分母に現れる素点は有限個（`Support.finite`）で、`v_p(y) ≤ 1` だからいずれも `p` と異なる。
中国剰余（`exists_forall_sub_mem_ideal`）で「`p` では `1`、分母の素点では十分高い冪」に
合同な `s` を取れば、`s·y` は全付値が `≤ 1` になるので整であり、
`y − s·y = (1 − s)·y` の `p`-付値は `< 1` である。 -/
theorem exists_integer_congr {L : Type} [Field L] [NumberField L]
    (p : HeightOneSpectrum (𝓞 L)) (y : L)
    (hy : HeightOneSpectrum.valuation L p y ≤ 1) :
    ∃ a : 𝓞 L, HeightOneSpectrum.valuation L p (y - algebraMap (𝓞 L) L a) < 1 := by
  rcases eq_or_ne y 0 with rfl | hy0
  · exact ⟨0, by simp⟩
  -- ★段 1: 分母に現れる素点は有限個
  have hfin := HeightOneSpectrum.Support.finite (𝓞 L) (K := L) y
  set T : Finset (HeightOneSpectrum (𝓞 L)) := hfin.toFinset with hT
  have hexp : ∀ q : HeightOneSpectrum (𝓞 L), ∃ m : ℕ,
      HeightOneSpectrum.valuation L q y ≤ WithZero.exp (m : ℤ) := by
    intro q
    have hne : HeightOneSpectrum.valuation L q y ≠ 0 := by simpa using hy0
    obtain ⟨k, hk⟩ : ∃ k : ℤ, HeightOneSpectrum.valuation L q y = WithZero.exp k :=
      ⟨WithZero.log _, (WithZero.exp_log hne).symm⟩
    exact ⟨k.toNat, by rw [hk, WithZero.exp_le_exp]; omega⟩
  choose n hn using hexp
  -- ★段 3: 中国剰余
  obtain ⟨s₀, hs₀⟩ := IsDedekindDomain.exists_forall_sub_mem_ideal
    (s := insert p T)
    (P := fun q : HeightOneSpectrum (𝓞 L) => q.asIdeal)
    (e := fun q => if q = p then 1 else n q)
    (fun q _ => HeightOneSpectrum.prime q)
    (fun q _ r _ hqr hc => hqr (HeightOneSpectrum.ext hc))
    (fun q => if (q : HeightOneSpectrum (𝓞 L)) = p then 1 else 0)
  have h1 : s₀ - 1 ∈ p.asIdeal ^ 1 := by
    have := hs₀ p (Finset.mem_insert_self p T)
    simpa using this
  -- ★段 4: `s₀·y` は整
  have hint : ∀ q : HeightOneSpectrum (𝓞 L),
      HeightOneSpectrum.valuation L q (algebraMap (𝓞 L) L s₀ * y) ≤ 1 := by
    intro q
    rw [map_mul]
    by_cases hq : q ∈ T
    · have hqp : q ≠ p := by
        intro hc
        rw [hc, hT, Set.Finite.mem_toFinset] at hq
        exact absurd hy (not_le.2 hq)
      have hmem : s₀ ∈ q.asIdeal ^ (n q) := by
        have hthis := hs₀ q (Finset.mem_insert_of_mem hq)
        simp only [hqp, if_false, sub_zero] at hthis
        exact hthis
      have hs : HeightOneSpectrum.valuation L q (algebraMap (𝓞 L) L s₀)
          ≤ WithZero.exp (-(n q : ℤ)) := by
        rw [HeightOneSpectrum.valuation_of_algebraMap]
        exact (HeightOneSpectrum.intValuation_le_pow_iff_mem q s₀ (n q)).mpr hmem
      calc HeightOneSpectrum.valuation L q (algebraMap (𝓞 L) L s₀)
              * HeightOneSpectrum.valuation L q y
          ≤ WithZero.exp (-(n q : ℤ)) * WithZero.exp ((n q : ℤ)) := mul_le_mul' hs (hn q)
        _ = 1 := by rw [← WithZero.exp_add]; simp
    · have hy1 : HeightOneSpectrum.valuation L q y ≤ 1 := by
        rw [hT, Set.Finite.mem_toFinset] at hq
        exact not_lt.1 hq
      have hs1 : HeightOneSpectrum.valuation L q (algebraMap (𝓞 L) L s₀) ≤ 1 :=
        HeightOneSpectrum.valuation_le_one q s₀
      calc HeightOneSpectrum.valuation L q (algebraMap (𝓞 L) L s₀)
              * HeightOneSpectrum.valuation L q y ≤ 1 * 1 := mul_le_mul' hs1 hy1
        _ = 1 := one_mul 1
  obtain ⟨a, ha⟩ := HeightOneSpectrum.mem_integers_of_valuation_le_one L _ hint
  refine ⟨a, ?_⟩
  rw [ha]
  have hfac : y - algebraMap (𝓞 L) L s₀ * y = (algebraMap (𝓞 L) L (1 - s₀)) * y := by
    rw [map_sub, map_one]; ring
  rw [hfac, map_mul]
  have hp1 : HeightOneSpectrum.valuation L p (algebraMap (𝓞 L) L (1 - s₀)) < 1 := by
    rw [HeightOneSpectrum.valuation_of_algebraMap,
      HeightOneSpectrum.intValuation_lt_one_iff_dvd]
    refine Ideal.dvd_span_singleton.2 ?_
    have hmem : s₀ - 1 ∈ p.asIdeal := by simpa using h1
    have hneg := neg_mem hmem
    simpa using hneg
  calc HeightOneSpectrum.valuation L p (algebraMap (𝓞 L) L (1 - s₀))
          * HeightOneSpectrum.valuation L p y
      ≤ HeightOneSpectrum.valuation L p (algebraMap (𝓞 L) L (1 - s₀)) * 1 :=
        mul_le_mul' le_rfl hy
    _ < 1 := by rw [mul_one]; exact hp1

/-! ## ★出典の紐付け(`.src`) -/

def exists_integer_congr.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(v_p(y) ≤ 1 なら y は 𝓞 L の元と p を法として合同。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GaloisRep
