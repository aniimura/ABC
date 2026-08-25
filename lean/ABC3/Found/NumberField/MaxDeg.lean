/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.NumberField.SplitSeparable
import ABC3.Found.NumberField.SplitExponent
import ABC3.Found.Divisor.ArithDivisor

/-!
# `[L:ℚ] = max_v deg(L,v)`(鎖 `cheb` の `cheb-max-deg`)

原典: S. Mochizuki, *The Geometry of Frobenioids I* [FrdI]、物理 p.116。

原文 (FrdI p.116):
> for the number of valuations ∈V(Li), including vi, that lie over the same valuation

原文 (FrdI p.116):
> of Q as vi. Then by Tchebotarev’s density theorem [cf., e.g., [Lang2], Chapter VIII,

## ★★原文の `deg(L, v)`

原文は `deg(L_i, v_i)` を「`v_i` を含め、`v_i` と**同じ `ℚ` の付値の上にある**
`V(L_i)` の付値の個数」と定める。`V(L) = FinitePlace L ⊕ InfinitePlace L`
(`ArithPlace`、`ArithDivisor.lean`)なので、これは 2 通りである:

* `v` が非アルキメデス(素数 `p` の上) —— `#{𝔭 ∣ p}` = `degAtPrime L p`
* `v` がアルキメデス —— `#(InfinitePlace L)` = `degAtInfinity L`

## ★★★Tchebotarev は要らない

原文は「Tchebotarev の密度定理により `[L:ℚ]` は `deg(L,v)` の最大値」と書く。
★**上界**は基本等式 `Σ_𝔭 e_𝔭 f_𝔭 = [L:ℚ]`(非アルキメデス)と
`r₁ + 2r₂ = [L:ℚ]`(アルキメデス)だけで出る。
★★**達成**は「完全分解する素数が 1 つある」ことだけで足り、それは
`SplitSeparable.lean` が **Chebotarev を使わずに**与えている。

## ★本ファイルで閉じること

| 定理 | 中身 |
|---|---|
| `degAtPrime_le` | `#{𝔭 ∣ p} ≤ [L:ℚ]` |
| `degAtInfinity_le` | `#(InfinitePlace L) ≤ [L:ℚ]` |
| `exists_degAtPrime_eq_finrank` | ★完全分解する素数で最大値が達成される |
| `finrank_isGreatest_deg` | ★★★★**`[L:ℚ] = max_v deg(L,v)`** |
-/

namespace ABC3.Found.NF

open _root_.NumberField Ideal
open scoped _root_.NumberField

variable (L : Type*) [Field L] [NumberField L]

/-! ## ★1. `deg(L, v)` の 2 つの形 -/

/-- ★**非アルキメデスの `deg(L,v)`** —— `p` の上の素イデアルの個数。 -/
noncomputable def degAtPrime (p : ℕ) : ℕ :=
  (primesOver (span {(p : ℤ)}) (𝓞 L)).ncard

/-- ★**アルキメデスの `deg(L,v)`** —— 無限素点の個数
(すべて `ℚ` の唯一のアルキメデス付値の上にある)。 -/
noncomputable def degAtInfinity : ℕ := Fintype.card (InfinitePlace L)

/-! ## ★2. 上界 -/

variable {L}

/-- ★`e_𝔭 · f_𝔭 ≥ 1`(底が `ℤ` の形)。 -/
theorem one_le_ramificationIdx_mul_inertiaDeg_int {p : ℕ} [Fact p.Prime]
    (P : Ideal (𝓞 L))
    (hP : P ∈ IsDedekindDomain.primesOverFinset (span {(p : ℤ)}) (𝓞 L)) :
    1 ≤ ramificationIdx (span {(p : ℤ)}) P * inertiaDeg (span {(p : ℤ)}) P := by
  have hne : span {(p : ℤ)} ≠ (⊥ : Ideal ℤ) := by simp [NeZero.ne p]
  rw [← Finset.mem_coe, IsDedekindDomain.coe_primesOverFinset hne (𝓞 L)] at hP
  haveI : P.IsPrime := hP.1
  haveI : P.LiesOver (span {(p : ℤ)}) := hP.2
  haveI : (span {(p : ℤ)}).IsMaximal := Int.ideal_span_isMaximal_of_prime p
  have h2 : 0 < inertiaDeg (span {(p : ℤ)}) P := Ideal.inertiaDeg_pos' _ P
  have h1 : ramificationIdx (span {(p : ℤ)}) P ≠ 0 :=
    Ideal.IsDedekindDomain.ramificationIdx_ne_zero_of_liesOver P hne
  exact Nat.one_le_iff_ne_zero.mpr (Nat.mul_ne_zero h1 h2.ne')

/-- ★★**非アルキメデスの上界** —— 基本等式 `Σ_𝔭 e_𝔭 f_𝔭 = [L:ℚ]` から。 -/
theorem degAtPrime_le (p : ℕ) [Fact p.Prime] :
    degAtPrime L p ≤ Module.finrank ℚ L := by
  classical
  have hne : span {(p : ℤ)} ≠ (⊥ : Ideal ℤ) := by simp [NeZero.ne p]
  have hsum := Ideal.sum_ramification_inertia (𝓞 L) ℚ L hne
  rw [degAtPrime, ← IsDedekindDomain.coe_primesOverFinset hne (𝓞 L), Set.ncard_coe_finset,
    ← hsum]
  have h := Finset.card_nsmul_le_sum
    (IsDedekindDomain.primesOverFinset (span {(p : ℤ)}) (𝓞 L))
    (fun P => ramificationIdx (span {(p : ℤ)}) P * inertiaDeg (span {(p : ℤ)}) P) 1
    (fun P hP => one_le_ramificationIdx_mul_inertiaDeg_int P hP)
  simpa using h

/-- ★★**アルキメデスの上界** —— `r₁ + 2r₂ = [L:ℚ]` から。 -/
theorem degAtInfinity_le : degAtInfinity L ≤ Module.finrank ℚ L := by
  have h1 := InfinitePlace.card_eq_nrRealPlaces_add_nrComplexPlaces L
  have h2 := InfinitePlace.card_add_two_mul_card_eq_rank L
  rw [degAtInfinity, h1]
  omega

/-! ## ★3. 最大値の達成 -/

/-- ★★★**完全分解する素数で最大値が達成される**。

★`SplitSeparable.lean` の `infinite_splitsCompletely'`(Chebotarev を使わない)から。 -/
theorem exists_degAtPrime_eq_finrank [IsGalois ℚ L] (θ : 𝓞 L)
    (hexp : _root_.RingOfIntegers.exponent θ ≠ 0) :
    ∃ p : ℕ, p.Prime ∧ degAtPrime L p = Module.finrank ℚ L := by
  obtain ⟨p, hp1, hp2⟩ := (infinite_splitsCompletely' θ hexp).nonempty
  exact ⟨p, hp1, hp2⟩

/-! ## ★4. 主定理 -/

/-- ★`deg(L, v)` が取りうる値の集合。 -/
def degSet (L : Type*) [Field L] [NumberField L] : Set ℕ :=
  {n | (∃ p : ℕ, p.Prime ∧ degAtPrime L p = n) ∨ n = degAtInfinity L}

/-- ★★★★★★**[FrdI] Theorem 6.4, (iv) の (a)** —— `[L:ℚ]` は `deg(L,v)` の**最大値**。

★★**Tchebotarev の密度定理は使わない** ——
上界は基本等式、達成は `infinite_splitsCompletely'` である。

原文 (FrdI p.116):
> for the number of valuations ∈V(Li), including vi, that lie over the same valuation -/
theorem finrank_isGreatest_deg [IsGalois ℚ L] (θ : 𝓞 L)
    (hexp : _root_.RingOfIntegers.exponent θ ≠ 0) :
    IsGreatest (degSet L) (Module.finrank ℚ L) := by
  obtain ⟨p, hp, hpe⟩ := exists_degAtPrime_eq_finrank θ hexp
  refine ⟨Or.inl ⟨p, hp, hpe⟩, ?_⟩
  rintro n (⟨q, hq, rfl⟩ | rfl)
  · haveI : Fact q.Prime := ⟨hq⟩
    exact degAtPrime_le q
  · exact degAtInfinity_le


/-! ## ★★★★★仮定を落とした形（2026-08-25）

`SplitExponent.lean` の `exists_exponent_ne_zero`（判別式で `hexp` を回収）により、
★以下は **`θ` も `hexp` も要らない**。 -/

variable (L) in
/-- ★★★**完全分解する素数で最大値が達成される**（★仮定なし）。 -/
theorem exists_degAtPrime_eq_finrank' [IsGalois ℚ L] :
    ∃ p : ℕ, p.Prime ∧ degAtPrime L p = Module.finrank ℚ L := by
  obtain ⟨θ, hexp⟩ := exists_exponent_ne_zero L
  exact exists_degAtPrime_eq_finrank θ hexp

variable (L) in
/-- ★★★★★★**[FrdI] Theorem 6.4, (iv) の (a)** —— `[L:ℚ]` は `deg(L,v)` の**最大値**
（★仮定なし）。

★★**Tchebotarev の密度定理は 1 度も使わない**。 -/
theorem finrank_isGreatest_deg' [IsGalois ℚ L] :
    IsGreatest (degSet L) (Module.finrank ℚ L) := by
  obtain ⟨θ, hexp⟩ := exists_exponent_ne_zero L
  exact finrank_isGreatest_deg θ hexp

/-! ### ★出典の紐付け -/

/-- ★★★★locator —— `Theorem 6.4, (iv)` の (a)(`[L:ℚ] = max_v deg(L,v)`)。 -/
def finrank_isGreatest_deg.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — [L:ℚ] は deg(L,v) の最大値(Tchebotarev を使わない)",
    sectionId := "frdi-thm-6-4" }


/-- ★★★★locator —— `Theorem 6.4, (iv)` の (a)(仮定なし)。 -/
def finrank_isGreatest_deg'.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 116,
    item := "Theorem 6.4, (iv) — [L:ℚ] は deg(L,v) の最大値(仮定なし・Tchebotarev 不使用)",
    sectionId := "frdi-thm-6-4" }

def finrank_isGreatest_deg'.needs : List ABC3.Meta.ProofObligation :=
  [ .citation "[ABC3]" "exists_exponent_ne_zero(判別式で hexp を回収)"
      (.inProject "ABC3" "ABC3.Found.NF.exists_exponent_ne_zero") 116,
    .citation "[ABC3]" "finrank_isGreatest_deg(θ・hexp つきの形)"
      (.inProject "ABC3" "ABC3.Found.NF.finrank_isGreatest_deg") 116 ]

end ABC3.Found.NF
