import ABC3.Found.NCBelyi.Normalize
import ABC3.Found.NCBelyi.Lemma21
import ABC3.Found.NCBelyi.BelyiComp

/-!
# [NCBelyi] Lemma 2.2 —— **基底段**(`Found`)

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.3–p.4。
**260 dpi 目視確認 2026-08-17**。

原文 (NCBelyi p.4):
> λ. Then, so long as |S| ≥ 4, the polynomial “f(x) + f0” of Lemma 2.1 determines

## ★★★何を取るか

`Lemma 2.2` は `|S|` についての帰納法である。★原文は
『**so long as `|S| ≥ 4`**』と書くだけで、**基底段(`|S| ≤ 3`)を書いていない**。

★★`S` は `0`, `∞` を含むので `|S| ≤ 3` ⟺ `|S\{0,∞}| ≤ 1` であり、
そこでは **1 次式 `x ↦ c·x`** で足りる(`Normalize.lean` の `exists_base_scale` が `c` を与える)。

★★★**1 次式は臨界点を持たない**(`(cx)′ = c ≠ 0`)ので、条件 (c) は**空虚に成り立つ**。
これが基底段が易しい理由である。

## ★結論の形

`ℙ¹` の `∞` は多項式なら常に `∞` に行くので、原文の (a)(b)(c) は
アフィン部分だけで書ける:

| 原文 | ここでの形 |
|---|---|
| (a) `φ(S) ⊆ {0,1,∞}` | `∀ α ∈ S, f(α) ∈ {0,1}` |
| (b) `φ(β) ∉ {0,1,∞}` | `f(β) ≠ 0` かつ `f(β) ≠ 1` |
| (c) `ℙ¹\{0,1,∞}` の上で不分岐 | 臨界**値**が `{0,1}` に入る |
-/

namespace ABC3.Found.NCBelyi

open Polynomial

/-! ## ★★基底段 -/

/-- ★★★★**[NCBelyi] Lemma 2.2 の基底段** —— `|S\{0}| ≤ 1` なら 1 次式で足りる。

原文 (NCBelyi p.4):
> λ. Then, so long as |S| ≥ 4, the polynomial “f(x) + f0” of Lemma 2.1 determines

★原文が書いていない段である。★★`Normalize.lean` の `exists_base_scale` が
係数 `c` を与え、`f = c·X` が (a)(b)(c) をすべて満たす。 -/
theorem lemma_2_2_base (S : Finset ℚ) (β : ℚ)
    (hpos : ∀ α ∈ S, α ≠ 0 → 0 < α)
    (hβ : β ∉ S) (hβ0 : β ≠ 0)
    (hcard : (S.erase 0).card ≤ 1) :
    ∃ f : ℚ[X], 0 < f.natDegree
      ∧ (∀ α ∈ S, f.eval α = 0 ∨ f.eval α = 1)
      ∧ f.eval β ≠ 0 ∧ f.eval β ≠ 1
      ∧ (∀ x : ℂ, (derivative (f.map (algebraMap ℚ ℂ))).eval x = 0 →
          (f.map (algebraMap ℚ ℂ)).eval x = 0 ∨ (f.map (algebraMap ℚ ℂ)).eval x = 1) := by
  classical
  have hposT : ∀ α ∈ S.erase 0, 0 < α := by
    intro α hα
    exact hpos α (Finset.mem_of_mem_erase hα) (Finset.ne_of_mem_erase hα)
  have hβT : β ∉ S.erase 0 := fun h => hβ (Finset.mem_of_mem_erase h)
  obtain ⟨c, hc0, hcT, hcβ0, hcβ1⟩ := exists_base_scale (S.erase 0) hcard hposT β hβ0 hβT
  refine ⟨C c * X, ?_, ?_, ?_, ?_, ?_⟩
  · rw [natDegree_C_mul hc0, natDegree_X]
    exact Nat.one_pos
  · intro α hα
    by_cases h : α = 0
    · left; simp [h]
    · right
      have : α ∈ S.erase 0 := Finset.mem_erase.2 ⟨h, hα⟩
      simpa using hcT α this
  · simpa using hcβ0
  · simpa using hcβ1
  · intro x hx
    exfalso
    rw [Polynomial.map_mul, Polynomial.map_C, Polynomial.map_X, derivative_mul,
      derivative_C, derivative_X] at hx
    simp only [zero_mul, mul_one, zero_add, eval_C] at hx
    exact hc0 ((map_eq_zero_iff (algebraMap ℚ ℂ) (algebraMap ℚ ℂ).injective).mp hx)

/-! ## ★★★★帰納段の部品 —— スケーリングとの合成

原文 (NCBelyi p.4):
> λ. Then, so long as |S| ≥ 4, the polynomial “f(x) + f0” of Lemma 2.1 determines

★帰納段は `Lemma 2.1` を **`λ·S`** に適用する。多項式の側では
`g(λ·x)`、すなわち **`g.comp (C λ * X)`** を扱うことになる。
★★その 2 つの基本性質(評価と臨界点)をここで取る。
-/

/-- ★**スケーリングとの合成の評価** —— `(g ∘ (λ·))(x) = g(λ x)`。 -/
theorem eval_comp_scale (g : ℚ[X]) (lam x : ℚ) :
    (g.comp (C lam * X)).eval x = g.eval (lam * x) := by
  rw [eval_comp]
  simp

/-- ★★**スケーリングは臨界点を `1/λ` 倍に移すだけ** ——
`(g(λ x))′ = λ·g′(λ x)` で `λ ≠ 0` だから。

★★★これがあれば、`Lemma 2.1` が与える「臨界点は `{0, 1, r}`」が
そのまま `λ·S` の側へ移る。 -/
theorem derivative_comp_scale_eval (g : ℚ[X]) (lam : ℚ) (x : ℂ) :
    (derivative ((g.comp (C lam * X)).map (algebraMap ℚ ℂ))).eval x
      = algebraMap ℚ ℂ lam * (derivative (g.map (algebraMap ℚ ℂ))).eval
          (algebraMap ℚ ℂ lam * x) := by
  rw [Polynomial.map_comp, derivative_comp]
  rw [Polynomial.map_mul, Polynomial.map_C, Polynomial.map_X, eval_mul, eval_comp]
  rw [derivative_mul, derivative_C, derivative_X]
  simp only [zero_mul, mul_one, zero_add, eval_C, eval_mul, eval_X]

/-- ★★★**臨界点の対応**(`λ ≠ 0`)。 -/
theorem derivative_comp_scale_eq_zero_iff (g : ℚ[X]) (lam : ℚ) (hlam : lam ≠ 0) (x : ℂ) :
    (derivative ((g.comp (C lam * X)).map (algebraMap ℚ ℂ))).eval x = 0
      ↔ (derivative (g.map (algebraMap ℚ ℂ))).eval (algebraMap ℚ ℂ lam * x) = 0 := by
  rw [derivative_comp_scale_eval]
  constructor
  · intro h
    rcases mul_eq_zero.mp h with h1 | h2
    · exact absurd ((map_eq_zero_iff (algebraMap ℚ ℂ)
        (algebraMap ℚ ℂ).injective).mp h1) hlam
    · exact h2
  · intro h
    rw [h, mul_zero]

/-- ★**スケーリングとの合成の次数**は元の次数に等しい(`λ ≠ 0`)。 -/
theorem natDegree_comp_scale (g : ℚ[X]) (lam : ℚ) (hlam : lam ≠ 0) :
    (g.comp (C lam * X)).natDegree = g.natDegree := by
  rw [natDegree_comp, natDegree_C_mul hlam, natDegree_X, mul_one]

/-! ## ★★★`Lemma 2.1` を `ℚ` 側へ移すための道具

`Lemma 2.1`(`Lemma21.lean`)は `ℝ` の上で書いてあるが、
`Lemma 2.2` の帰納法は `ℚ` の上で回る。
★差は `f₀ ≙ -inf'` のところだけである——埋め込みが単調なので下確は保たれる。
-/

/-- ★★**`ℚ → ℝ` の埋め込みは `inf'` と可換**である。

★単調な埋め込みなので両向きの不等式が出る。
★★`Lemma 2.1` の `f₀ ≙ -S.inf'` を `ℚ` 側で定めても同じものになる、ということである。 -/
theorem rat_cast_inf' {ι : Type*} (S : Finset ι) (hne : S.Nonempty) (g : ι → ℚ) :
    ((S.inf' hne g : ℚ) : ℝ) = S.inf' hne (fun i => (g i : ℝ)) := by
  refine le_antisymm ?_ ?_
  · refine Finset.le_inf' hne _ (fun i hi => ?_)
    exact_mod_cast Finset.inf'_le g hi
  · obtain ⟨j, hjS, hj⟩ := Finset.exists_mem_eq_inf' hne g
    rw [hj]
    exact Finset.inf'_le (fun i => (g i : ℝ)) hjS

/-! ## ★★★★`ℚ` 側の `f` と `f₀` -/

/-- ★`Lemma 2.1` の `f(x) = x^{a+1}(x−1)^{b+1}` を `ℚ` で書いたもの。 -/
def belyiVal (a b : ℕ) (x : ℚ) : ℚ := x ^ (a + 1) * (x - 1) ^ (b + 1)

/-- ★`Lemma 2.1` の `f₀ ≙ -inf'` を `ℚ` で書いたもの。 -/
noncomputable def belyiShift (a b : ℕ) (S : Finset ℚ) (hne : S.Nonempty) : ℚ :=
  - S.inf' hne (belyiVal a b)

/-- ★★**`ℚ` で定めた `f` をキャストすると `ℝ` 側の `f` になる**。 -/
theorem cast_belyiVal (a b : ℕ) (x : ℚ) :
    ((belyiVal a b x : ℚ) : ℝ) = (x : ℝ) ^ (a + 1) * ((x : ℝ) - 1) ^ (b + 1) := by
  simp [belyiVal]

/-- ★★★**`ℚ` で定めた `f₀` をキャストすると `ℝ` 側の `f₀` になる**。

★像の上の `inf'` を `Finset.inf'_image` で元の集合の `inf'` に戻し、
`rat_cast_inf'` でキャストを外に出す。 -/
theorem cast_inf'_image (a b : ℕ) (S : Finset ℚ) (hne : S.Nonempty) :
    ((S.image (fun q : ℚ => (q : ℝ))).inf' (hne.image _)
        (fun x : ℝ => x ^ (a + 1) * (x - 1) ^ (b + 1)))
      = ((S.inf' hne (belyiVal a b) : ℚ) : ℝ) := by
  rw [Finset.inf'_image, rat_cast_inf' S hne (belyiVal a b)]
  refine Finset.inf'_congr hne rfl (fun x _ => ?_)
  simp [Function.comp, belyiVal]

/-! ## ★★★★★`Lemma 2.1` の (c)(d) を `ℚ` へ戻す -/

/-- ★★★★★**`Lemma 2.1` の (c)(d) の `ℚ` 版**。

原文 (NCBelyi p.2):
> Lemma 2.1. (Separating Properties of Belyi Maps) Let C ∈ R be such

★`Lemma21.lean` の `lemma_2_1` は `ℝ` の上で書いてあるが、
`Lemma 2.2` の帰納法は `ℚ` の上で回る。埋め込みが単射・単調なので
結論はそのまま戻る。

* (c) `f(β) ∉ f(S)`
* (d) `(f(β)+f₀)/(f(α)+f₀) ≥ 2` —— ★これが帰納仮定の条件 (iii) を再現する -/
theorem lemma_2_1_rat (a b : ℕ) (S : Finset ℚ) (hne : S.Nonempty) (r β : ℚ)
    (hr : (r : ℝ) = ((a : ℝ) + 1) / (((a : ℝ) + 1) + ((b : ℝ) + 1)))
    (h0S : (0 : ℚ) ∈ S) (h1S : (1 : ℚ) ∈ S) (hrS : r ∈ S)
    (hS : ∀ α ∈ S, α = 0 ∨ α = r ∨ α = 1 ∨ 1 < α)
    (hβ : ∀ α ∈ S, α ≠ 0 → 2 * α ≤ β) :
    (∀ α ∈ S, belyiVal a b β ≠ belyiVal a b α)
  ∧ (∀ α ∈ S, belyiVal a b α + belyiShift a b S hne ≠ 0 →
      2 ≤ (belyiVal a b β + belyiShift a b S hne)
          / (belyiVal a b α + belyiShift a b S hne)) := by
  classical
  have hne' : (S.image (fun q : ℚ => (q : ℝ))).Nonempty := hne.image _
  have hmem : ∀ α ∈ S, ((α : ℝ)) ∈ (S.image (fun q : ℚ => (q : ℝ))) := fun α hα => Finset.mem_image.2 ⟨α, hα, rfl⟩
  have h0img : (0 : ℝ) ∈ (S.image (fun q : ℚ => (q : ℝ))) := by simpa using hmem 0 h0S
  have h1img : (1 : ℝ) ∈ (S.image (fun q : ℚ => (q : ℝ))) := by simpa using hmem 1 h1S
  have hrimg : (r : ℝ) ∈ (S.image (fun q : ℚ => (q : ℝ))) := hmem r hrS
  have hS'' : ∀ α ∈ (S.image (fun q : ℚ => (q : ℝ))), α = 0 ∨ α = (r : ℝ) ∨ α = 1 ∨ 1 < α := by
    intro α hα
    obtain ⟨q, hq, rfl⟩ := Finset.mem_image.1 hα
    rcases hS q hq with h | h | h | h
    · left; exact_mod_cast congrArg (fun z : ℚ => (z : ℝ)) h
    · right; left; exact_mod_cast congrArg (fun z : ℚ => (z : ℝ)) h
    · right; right; left; exact_mod_cast congrArg (fun z : ℚ => (z : ℝ)) h
    · right; right; right; exact_mod_cast h
  have hβ' : ∀ α ∈ (S.image (fun q : ℚ => (q : ℝ))), α ≠ 0 → (2 : ℝ) * α ≤ (β : ℝ) := by
    intro α hα hne0
    obtain ⟨q, hq, rfl⟩ := Finset.mem_image.1 hα
    have hq0 : q ≠ 0 := by
      intro h
      exact hne0 (by rw [h]; norm_num)
    exact_mod_cast hβ q hq hq0
  obtain ⟨-, -, hc, hd⟩ :=
    lemma_2_1 2 (le_refl 2) a b (S.image (fun q : ℚ => (q : ℝ))) hne' (r : ℝ) (β : ℝ) hr h0img h1img hrimg hS'' hβ'
  have hinf := cast_inf'_image a b S hne
  refine ⟨?_, ?_⟩
  · intro α hα hcon
    refine hc ((α : ℝ)) (hmem α hα) ?_
    rw [← cast_belyiVal, ← cast_belyiVal, hcon]
  · intro α hα hnz
    have hnz' : ((α : ℝ)) ^ (a + 1) * (((α : ℝ)) - 1) ^ (b + 1)
        - (S.image (fun q : ℚ => (q : ℝ))).inf' hne' (fun x : ℝ => x ^ (a + 1) * (x - 1) ^ (b + 1)) ≠ 0 := by
      rw [← cast_belyiVal, hinf]
      intro h
      refine hnz ?_
      have : ((belyiVal a b α + belyiShift a b S hne : ℚ) : ℝ) = 0 := by
        rw [belyiShift]
        push_cast
        linarith [h]
      exact_mod_cast this
    have hdd := hd ((α : ℝ)) (hmem α hα) hnz'
    rw [← cast_belyiVal, ← cast_belyiVal, hinf] at hdd
    have hcast : ((belyiVal a b β + belyiShift a b S hne : ℚ) : ℝ)
        = (belyiVal a b β : ℝ) - ((S.inf' hne (belyiVal a b) : ℚ) : ℝ) := by
      rw [belyiShift]; push_cast; ring
    have hcast2 : ((belyiVal a b α + belyiShift a b S hne : ℚ) : ℝ)
        = (belyiVal a b α : ℝ) - ((S.inf' hne (belyiVal a b) : ℚ) : ℝ) := by
      rw [belyiShift]; push_cast; ring
    have : (2 : ℝ) ≤ ((belyiVal a b β + belyiShift a b S hne : ℚ) : ℝ)
        / ((belyiVal a b α + belyiShift a b S hne : ℚ) : ℝ) := by
      rw [hcast, hcast2]; exact hdd
    exact_mod_cast this

/-! ## ★★★★★★帰納段の帳尻 —— `S₂` の 3 つの性質

原文 (NCBelyi p.4):
> λ. Then, so long as |S| ≥ 4, the polynomial “f(x) + f0” of Lemma 2.1 determines

帰納段は `S₂ ≙ (f + f₀)(λ·S)` を作って帰納仮定を適用する。
★`S₂` が帰納仮定の形を満たすことを、ここで 3 つに分けて取る。
-/

/-- ★**`f(0) = f(1) = 0`** —— `Lemma 2.1` の (a)。 -/
theorem belyiVal_zero (a b : ℕ) : belyiVal a b 0 = 0 := by
  simp [belyiVal]

theorem belyiVal_one (a b : ℕ) : belyiVal a b 1 = 0 := by
  simp [belyiVal]

/-- ★★★**`|S₂| < |S₁|`** —— `0` と `1` が同じ値へ行くから。

原文 (NCBelyi p.4):
> which the cardinalities of S , S satisfy |S | < |S|.

★原文は `|S′| < |S|` とだけ書くが、**その根拠は `f(0) = f(1)`** である。 -/
theorem card_image_belyiVal_lt (a b : ℕ) (S₁ : Finset ℚ) (f0 : ℚ)
    (h0 : (0 : ℚ) ∈ S₁) (h1 : (1 : ℚ) ∈ S₁) :
    (S₁.image (fun α => belyiVal a b α + f0)).card < S₁.card := by
  refine card_image_lt_of_collision h0 h1 (by norm_num) ?_
  rw [belyiVal_zero, belyiVal_one]

/-- ★★**`f₀` を足すと値は非負になる** —— `f₀ ≙ -inf'` だから。 -/
theorem belyiVal_add_shift_nonneg (a b : ℕ) (S₁ : Finset ℚ) (hne : S₁.Nonempty)
    {α : ℚ} (hα : α ∈ S₁) :
    0 ≤ belyiVal a b α + belyiShift a b S₁ hne := by
  rw [belyiShift]
  have := Finset.inf'_le (belyiVal a b) hα
  linarith

/-- ★★★**`0 ∈ S₂`** —— 下確を達成する点で `f + f₀` が `0` になる。 -/
theorem zero_mem_image_belyiVal_add_shift (a b : ℕ) (S₁ : Finset ℚ) (hne : S₁.Nonempty) :
    (0 : ℚ) ∈ S₁.image (fun α => belyiVal a b α + belyiShift a b S₁ hne) := by
  obtain ⟨j, hjS, hj⟩ := Finset.exists_mem_eq_inf' hne (belyiVal a b)
  refine Finset.mem_image.2 ⟨j, hjS, ?_⟩
  rw [belyiShift, ← hj]
  ring

/-! ## ★★★★★★★帰納段の多項式とその臨界点 -/

/-- ★帰納段で使う多項式 —— `Lemma 2.1` の `f + f₀` を `λ` 倍したもの。

原文 (NCBelyi p.4):
> λ. Then, so long as |S| ≥ 4, the polynomial “f(x) + f0” of Lemma 2.1 determines -/
noncomputable def stepPoly (a b : ℕ) (f0 lam : ℚ) : ℚ[X] :=
  ((X : ℚ[X]) ^ (a + 1) * (X - 1) ^ (b + 1) + C f0).comp (C lam * X)

/-- ★★**`stepPoly` の評価** —— `λ x` での `f` の値に `f₀` を足したもの。 -/
theorem eval_stepPoly (a b : ℕ) (f0 lam x : ℚ) :
    (stepPoly a b f0 lam).eval x = belyiVal a b (lam * x) + f0 := by
  rw [stepPoly, eval_comp]
  simp [belyiVal]

/-- ★**`stepPoly` は非定数**(`λ ≠ 0`)。 -/
theorem natDegree_stepPoly_pos (a b : ℕ) (f0 lam : ℚ) (hlam : lam ≠ 0) :
    0 < (stepPoly a b f0 lam).natDegree := by
  rw [stepPoly, natDegree_comp_scale _ _ hlam]
  have h1 : ((X : ℚ[X]) - 1).Monic := by simpa using monic_X_sub_C (1 : ℚ)
  have hdeg1 : ((X : ℚ[X]) - 1).natDegree = 1 := by simpa using natDegree_X_sub_C (1 : ℚ)
  have hdeg : ((X : ℚ[X]) ^ (a + 1) * (X - 1) ^ (b + 1)).natDegree = (a + 1) + (b + 1) := by
    rw [Monic.natDegree_mul (monic_X_pow (a + 1)) (h1.pow (b + 1)), natDegree_X_pow,
      natDegree_pow, hdeg1, mul_one]
  rw [natDegree_add_C, hdeg]
  omega

/-- ★★★★★**`stepPoly` の臨界点は `λ` 倍すると `{0, 1, r}` に入る**。

原文 (NCBelyi p.2):
> mediately from the fact that:

★`f₀` を足しても導関数は変わらず、スケーリングは臨界点を `1/λ` 倍に移すだけである
(第 399)。★★あとは `Lemma21.lean` の `belyi_critical` そのもの。 -/
theorem crit_stepPoly (a b : ℕ) (f0 lam : ℚ) (hlam : lam ≠ 0) (x : ℂ)
    (hx : (derivative ((stepPoly a b f0 lam).map (algebraMap ℚ ℂ))).eval x = 0) :
    algebraMap ℚ ℂ lam * x = 0 ∨ algebraMap ℚ ℂ lam * x = 1
      ∨ algebraMap ℚ ℂ lam * x = ((a : ℂ) + 1) / ((a : ℂ) + 1 + ((b : ℂ) + 1)) := by
  rw [stepPoly, derivative_comp_scale_eq_zero_iff _ _ hlam] at hx
  have hmap : (((X : ℚ[X]) ^ (a + 1) * (X - 1) ^ (b + 1) + C f0).map (algebraMap ℚ ℂ))
      = (X : ℂ[X]) ^ (a + 1) * (X - 1) ^ (b + 1) + C (algebraMap ℚ ℂ f0) := by
    simp [Polynomial.map_add, Polynomial.map_mul, Polynomial.map_pow, Polynomial.map_sub]
  rw [hmap, derivative_add, derivative_C, add_zero] at hx
  exact belyi_critical a b hx

/-! ## ★★★★★★★★帰納段の受け皿 —— `f + f₀` が条件 (i)(ii)(iii) を作り直す -/

/-- ★★★★★★**`Lemma 2.2` の帰納段の帳尻**。

原文 (NCBelyi p.4):
> some β, S that satisfy conditions (i), (ii), (iii) of the present Lemma 2.2, but for

★原文はこの 1 文で済ませるが、中身は
`Lemma 2.1` の (a)(c)(d) を条件 (i)(ii)(iii) へ**翻訳する**作業である。

| 出す性質 | 根拠 |
|---|---|
| `0 ∈ S′`(条件 (i)) | 下確を達成する点で `f + f₀ = 0` |
| `α ∈ S′\{0} → α > 0`(条件 (ii)) | `f₀ ≝ -inf'` だから `f + f₀ ≥ 0` |
| `β′ ∉ S′` | `Lemma 2.1` (c) |
| `β′/α ≥ 2`(条件 (iii)) | `Lemma 2.1` (d) |
| `|S′| < |S|` | `f(0) = f(1)` |
-/
theorem lemma_2_2_shift (a b : ℕ) (S₁ : Finset ℚ) (hne : S₁.Nonempty) (r β₁ : ℚ)
    (hr : (r : ℝ) = ((a : ℝ) + 1) / (((a : ℝ) + 1) + ((b : ℝ) + 1)))
    (h0S : (0 : ℚ) ∈ S₁) (h1S : (1 : ℚ) ∈ S₁) (hrS : r ∈ S₁)
    (hS : ∀ α ∈ S₁, α = 0 ∨ α = r ∨ α = 1 ∨ 1 < α)
    (hβ : ∀ α ∈ S₁, α ≠ 0 → 2 * α ≤ β₁) :
    (0 : ℚ) ∈ S₁.image (fun α => belyiVal a b α + belyiShift a b S₁ hne)
    ∧ (∀ α ∈ S₁.image (fun α => belyiVal a b α + belyiShift a b S₁ hne), α ≠ 0 → 0 < α)
    ∧ (belyiVal a b β₁ + belyiShift a b S₁ hne)
        ∉ S₁.image (fun α => belyiVal a b α + belyiShift a b S₁ hne)
    ∧ (∀ α ∈ S₁.image (fun α => belyiVal a b α + belyiShift a b S₁ hne), α ≠ 0 →
        2 * α ≤ belyiVal a b β₁ + belyiShift a b S₁ hne)
    ∧ (S₁.image (fun α => belyiVal a b α + belyiShift a b S₁ hne)).card < S₁.card := by
  classical
  obtain ⟨hc, hd⟩ := lemma_2_1_rat a b S₁ hne r β₁ hr h0S h1S hrS hS hβ
  refine ⟨zero_mem_image_belyiVal_add_shift a b S₁ hne, ?_, ?_, ?_, ?_⟩
  · intro α hα hα0
    obtain ⟨γ, hγ, rfl⟩ := Finset.mem_image.1 hα
    exact lt_of_le_of_ne (belyiVal_add_shift_nonneg a b S₁ hne hγ) (Ne.symm hα0)
  · intro hcon
    obtain ⟨γ, hγ, hγe⟩ := Finset.mem_image.1 hcon
    exact hc γ hγ (by linarith [hγe])
  · intro α hα hα0
    obtain ⟨γ, hγ, rfl⟩ := Finset.mem_image.1 hα
    have hpos : 0 < belyiVal a b γ + belyiShift a b S₁ hne :=
      lt_of_le_of_ne (belyiVal_add_shift_nonneg a b S₁ hne hγ) (Ne.symm hα0)
    have hdd := hd γ hγ hα0
    rwa [le_div_iff₀ hpos] at hdd
  · exact card_image_belyiVal_lt a b S₁ _ h0S h1S

/-! ## ★★★★★★★★★帰納段 -/

/-- ★`λ` 倍を打ち消して「臨界点は有理点である」を取り出す。 -/
theorem eq_ratCast_of_scale {lam s : ℚ} (hlam : lam ≠ 0) {x : ℂ}
    (h : algebraMap ℚ ℂ lam * x = algebraMap ℚ ℂ (lam * s)) : x = algebraMap ℚ ℂ s := by
  have hl : algebraMap ℚ ℂ lam ≠ 0 := fun hc =>
    hlam ((algebraMap ℚ ℂ).injective (by simpa using hc))
  rw [map_mul] at h
  exact mul_left_cancel₀ hl h

/-- ★★★★★★★★**[NCBelyi] Lemma 2.2 の帰納段**。

原文 (NCBelyi p.4):
> Indeed, we induct on the cardinality |S| of S and apply Lemma 2.1 [with,
> say, C = 2] to the set λ · S ⊆ P1

★★`λ` は `Normalize.lean` の `exists_normalizing_scale`(= `1/α₂`)、
`a, b` は `exists_num_den`(= `r` の分子と分母−分子)で決まる。
★★★これで `λ·S` が `Lemma 2.1` の仮定 `α ∈ {0, r, 1} ∪ (1,∞)` を**ちょうど**満たす。

出す多項式は `stepPoly a b f₀ λ = (f + f₀)(λx)` である。 -/
theorem lemma_2_2_step (S : Finset ℚ) (β : ℚ)
    (h0S : (0 : ℚ) ∈ S)
    (hpos : ∀ α ∈ S, α ≠ 0 → 0 < α)
    (hratio : ∀ α ∈ S, α ≠ 0 → 2 * α ≤ β)
    (hcard : 2 ≤ (S.erase 0).card) :
    ∃ (h : ℚ[X]) (S₂ : Finset ℚ) (β₂ : ℚ),
      0 < h.natDegree
      ∧ (0 : ℚ) ∈ S₂
      ∧ (∀ α ∈ S, h.eval α ∈ S₂)
      ∧ h.eval β = β₂
      ∧ (∀ α ∈ S₂, α ≠ 0 → 0 < α)
      ∧ β₂ ∉ S₂
      ∧ (∀ α ∈ S₂, α ≠ 0 → 2 * α ≤ β₂)
      ∧ (S₂.erase 0).card < (S.erase 0).card
      ∧ (∀ x : ℂ, (derivative (h.map (algebraMap ℚ ℂ))).eval x = 0 →
          ∃ s ∈ S, (h.map (algebraMap ℚ ℂ)).eval x = algebraMap ℚ ℂ (h.eval s)) := by
  classical
  have hposT : ∀ α ∈ S.erase 0, 0 < α := fun α hα =>
    hpos α (Finset.mem_of_mem_erase hα) (Finset.ne_of_mem_erase hα)
  obtain ⟨lam, r, hlam0, hr0, hr1, h1mem, hrmem, hall⟩ :=
    exists_normalizing_scale (S.erase 0) hposT hcard
  obtain ⟨a, b, hr⟩ := exists_num_den r hr0 hr1
  have hlamne : lam ≠ 0 := hlam0.ne'
  obtain ⟨S₁, hS₁def⟩ : ∃ S₁ : Finset ℚ, S₁ = S.image (fun t => lam * t) := ⟨_, rfl⟩
  have hmemS₁ : ∀ t ∈ S, lam * t ∈ S₁ := by
    intro t ht
    rw [hS₁def]
    exact Finset.mem_image.2 ⟨t, ht, rfl⟩
  have hsub : (S.erase 0).image (fun t => lam * t) ⊆ S₁ := by
    rw [hS₁def]
    exact Finset.image_subset_image (Finset.erase_subset _ _)
  have h0S₁ : (0 : ℚ) ∈ S₁ := by
    have := hmemS₁ 0 h0S
    rwa [mul_zero] at this
  have h1S₁ : (1 : ℚ) ∈ S₁ := hsub h1mem
  have hrS₁ : r ∈ S₁ := hsub hrmem
  have hne₁ : S₁.Nonempty := ⟨0, h0S₁⟩
  have hS₁all : ∀ α ∈ S₁, α = 0 ∨ α = r ∨ α = 1 ∨ 1 < α := by
    intro α hα
    rw [hS₁def] at hα
    obtain ⟨t, htS, rfl⟩ := Finset.mem_image.1 hα
    by_cases h : t = 0
    · left; rw [h, mul_zero]
    · exact Or.inr (hall _ (Finset.mem_image.2 ⟨t, Finset.mem_erase.2 ⟨h, htS⟩, rfl⟩))
  have hβ₁ : ∀ α ∈ S₁, α ≠ 0 → 2 * α ≤ lam * β := by
    intro α hα hα0
    rw [hS₁def] at hα
    obtain ⟨t, htS, rfl⟩ := Finset.mem_image.1 hα
    have ht0 : t ≠ 0 := by
      intro h
      exact hα0 (by rw [h, mul_zero])
    calc 2 * (lam * t) = lam * (2 * t) := by ring
      _ ≤ lam * β := mul_le_mul_of_nonneg_left (hratio t htS ht0) hlam0.le
  obtain ⟨h0S₂, hpos₂, hβ₂notmem, hratio₂, hcardlt⟩ :=
    lemma_2_2_shift a b S₁ hne₁ r (lam * β) hr h0S₁ h1S₁ hrS₁ hS₁all hβ₁
  refine ⟨stepPoly a b (belyiShift a b S₁ hne₁) lam,
    S₁.image (fun α => belyiVal a b α + belyiShift a b S₁ hne₁),
    belyiVal a b (lam * β) + belyiShift a b S₁ hne₁,
    natDegree_stepPoly_pos a b _ lam hlamne, h0S₂, ?_, eval_stepPoly a b _ lam β,
    hpos₂, hβ₂notmem, hratio₂, ?_, ?_⟩
  · intro α hα
    rw [eval_stepPoly]
    exact Finset.mem_image.2 ⟨lam * α, hmemS₁ α hα, rfl⟩
  · -- `|S₂\{0}| < |S\{0}|`
    have hcardS₁ : S₁.card = S.card := by
      rw [hS₁def, Finset.card_image_of_injective _ (mul_right_injective₀ hlamne)]
    have h1le : 1 ≤ (S₁.image (fun α => belyiVal a b α + belyiShift a b S₁ hne₁)).card :=
      Finset.card_pos.2 ⟨0, h0S₂⟩
    rw [Finset.card_erase_of_mem h0S₂, Finset.card_erase_of_mem h0S]
    omega
  · -- 臨界点は `λ·S` の中の `{0, 1, r}` から来る
    have hrC : algebraMap ℚ ℂ r = ((a : ℂ) + 1) / ((a : ℂ) + 1 + ((b : ℂ) + 1)) := by
      have hcast := congrArg (fun z : ℝ => (z : ℂ)) hr
      push_cast at hcast
      simpa using hcast
    intro x hx
    rcases crit_stepPoly a b (belyiShift a b S₁ hne₁) lam hlamne x hx with hcx | hcx | hcx
    · refine ⟨0, h0S, ?_⟩
      have hxs : x = algebraMap ℚ ℂ (0 : ℚ) := by
        refine eq_ratCast_of_scale hlamne ?_
        rw [hcx, mul_zero, map_zero]
      rw [hxs, eval_map_ratCast]
    · obtain ⟨t, htT, hts⟩ := Finset.mem_image.1 h1mem
      refine ⟨t, Finset.mem_of_mem_erase htT, ?_⟩
      have hxs : x = algebraMap ℚ ℂ t := by
        refine eq_ratCast_of_scale hlamne ?_
        rw [hcx, hts, map_one]
      rw [hxs, eval_map_ratCast]
    · obtain ⟨t, htT, hts⟩ := Finset.mem_image.1 hrmem
      refine ⟨t, Finset.mem_of_mem_erase htT, ?_⟩
      have hxs : x = algebraMap ℚ ℂ t := by
        refine eq_ratCast_of_scale hlamne ?_
        rw [hcx, hts, hrC]
      rw [hxs, eval_map_ratCast]

/-! ## ★★★★★★★★★★`Lemma 2.2` 本体 —— `|S\{0}|` についての帰納法 -/

/-- ★★★★★★★★★**帰納法の本体**(`n` は `|S\{0}|` の上界)。

原文 (NCBelyi p.4):
> hypothesis and composing the resulting morphisms P1

★★合成は `BelyiComp.lean` の `comp_crit_of_rel` が受ける
——★★★**合成される 2 つは対等でない**:
`stepPoly` は `{0,1}` を `{0,1}` へ写さない(`f₀` へ写す)ので `IsBelyiPoly` ではなく、
臨界値が `stepPoly(S)` に入るという**相対的な**性質しか持たない。
それを帰納法の仮定 `g(S₂) ⊆ {0,1}` が吸収する。 -/
theorem lemma_2_2_aux : ∀ (n : ℕ) (S : Finset ℚ) (β : ℚ),
    (S.erase 0).card ≤ n →
    (0 : ℚ) ∈ S →
    (∀ α ∈ S, α ≠ 0 → 0 < α) →
    β ∉ S → β ≠ 0 →
    (∀ α ∈ S, α ≠ 0 → 2 * α ≤ β) →
    ∃ f : ℚ[X], 0 < f.natDegree
      ∧ (∀ α ∈ S, f.eval α = 0 ∨ f.eval α = 1)
      ∧ f.eval β ≠ 0 ∧ f.eval β ≠ 1
      ∧ (∀ x : ℂ, (derivative (f.map (algebraMap ℚ ℂ))).eval x = 0 →
          (f.map (algebraMap ℚ ℂ)).eval x = 0 ∨ (f.map (algebraMap ℚ ℂ)).eval x = 1) := by
  intro n
  induction n with
  | zero =>
    intro S β hcard _ hpos hβ hβ0 _
    exact lemma_2_2_base S β hpos hβ hβ0 (by omega)
  | succ n ih =>
    intro S β hcard h0S hpos hβ hβ0 hratio
    by_cases hsmall : (S.erase 0).card ≤ 1
    · exact lemma_2_2_base S β hpos hβ hβ0 hsmall
    · push_neg at hsmall
      obtain ⟨h, S₂, β₂, hdeg, h0S₂, hmaps, hevalβ, hpos₂, hβ₂, hratio₂, hcard₂, hcrit⟩ :=
        lemma_2_2_step S β h0S hpos hratio (by omega)
      have hβ₂0 : β₂ ≠ 0 := by
        intro hc
        rw [hc] at hβ₂
        exact hβ₂ h0S₂
      obtain ⟨g, hgdeg, hgS, hgβ0, hgβ1, hgcrit⟩ :=
        ih S₂ β₂ (by omega) h0S₂ hpos₂ hβ₂ hβ₂0 hratio₂
      refine ⟨g.comp h, ?_, ?_, ?_, ?_, ?_⟩
      · rw [natDegree_comp]
        exact Nat.mul_pos hgdeg hdeg
      · intro α hα
        rw [eval_comp]
        exact hgS _ (hmaps α hα)
      · rw [eval_comp, hevalβ]
        exact hgβ0
      · rw [eval_comp, hevalβ]
        exact hgβ1
      · exact comp_crit_of_rel hcrit (fun s hs => hgS _ (hmaps s hs)) hgcrit

/-- ★★★★★★★★★★**[NCBelyi] Lemma 2.2**(Belyi Maps Noncritical at Prescribed Rational Points)。

原文 (NCBelyi p.3):
> Lemma 2.2.
> (Belyi Maps Noncritical at Prescribed Rational Points)

`S ⊆ ℙ¹(ℚ)` が (i) `0, ∞ ∈ S`、(ii) `α ∈ S\{0,∞} → α > 0` を満たし、
`β ∈ ℚ\S` が (iii) `β/α ≥ 2`(`∀ α ∈ S\{0,∞}`)を満たすとき、
**非定数多項式 `f ∈ ℚ[x]`** が存在して
(a) `φ(S) ⊆ {0,1,∞}`、(b) `φ(β) ∉ {0,1,∞}`、(c) `φ` は `ℙ¹\{0,1,∞}` 上不分岐。

★`∞` は多項式なら常に `∞` へ行くので、アフィン部分だけで書いてある
(`Lemma22.lean` 冒頭の対応表)。★条件 (iii) は `β/α ≥ 2` を
**割り算を避けて `2·α ≤ β`** と書いた(`α > 0` なので同値)。 -/
theorem lemma_2_2 (S : Finset ℚ) (β : ℚ)
    (h0S : (0 : ℚ) ∈ S)
    (hpos : ∀ α ∈ S, α ≠ 0 → 0 < α)
    (hβ : β ∉ S) (hβ0 : β ≠ 0)
    (hratio : ∀ α ∈ S, α ≠ 0 → 2 * α ≤ β) :
    ∃ f : ℚ[X], 0 < f.natDegree
      ∧ (∀ α ∈ S, f.eval α = 0 ∨ f.eval α = 1)
      ∧ f.eval β ≠ 0 ∧ f.eval β ≠ 1
      ∧ (∀ x : ℂ, (derivative (f.map (algebraMap ℚ ℂ))).eval x = 0 →
          (f.map (algebraMap ℚ ℂ)).eval x = 0 ∨ (f.map (algebraMap ℚ ℂ)).eval x = 1) :=
  lemma_2_2_aux (S.erase 0).card S β le_rfl h0S hpos hβ hβ0 hratio

/-! ## ★出典の紐付け(`.src`) -/

def lemma_2_2_base.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.2(基底段——原文が書いていない |S| ≤ 3 の場合)",
    sectionId := "ncbelyi-lemma-2-2" }

def derivative_comp_scale_eq_zero_iff.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.2(帰納段——スケーリングは臨界点を移すだけ)",
    sectionId := "ncbelyi-lemma-2-2" }

def rat_cast_inf'.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.2(帰納段——Lemma 2.1 を ℚ 側へ移すための inf' のキャスト)",
    sectionId := "ncbelyi-lemma-2-2" }

def belyiVal.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 2,
    item := "Lemma 2.1 の f(x) = x^{a+1}(x−1)^{b+1}(ℚ 版)",
    sectionId := "ncbelyi-lemma-2-1" }

def cast_inf'_image.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.2(帰納段——f₀ を ℚ で定めても同じものになる)",
    sectionId := "ncbelyi-lemma-2-2" }

def lemma_2_1_rat.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 2,
    item := "Lemma 2.1, (c)(d)(ℚ 版——Lemma 2.2 の帰納法が使う形)",
    sectionId := "ncbelyi-lemma-2-1" }

def card_image_belyiVal_lt.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.2(帰納段——f(0) = f(1) から |S′| < |S|)",
    sectionId := "ncbelyi-lemma-2-2" }

def zero_mem_image_belyiVal_add_shift.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.2(帰納段——f₀ を足すと下確が 0 になる)",
    sectionId := "ncbelyi-lemma-2-2" }

def stepPoly.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.2(帰納段の多項式——f + f₀ を λ 倍したもの)",
    sectionId := "ncbelyi-lemma-2-2" }

def crit_stepPoly.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 2,
    item := "Lemma 2.1, (b)(スケーリングと f₀ を通した形)",
    sectionId := "ncbelyi-lemma-2-1" }

def lemma_2_2_shift.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.2(帰納段——β′, S′ が条件 (i)(ii)(iii) を満たす)",
    sectionId := "ncbelyi-lemma-2-2" }

def lemma_2_2_step.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.2(帰納段——λ·S へ Lemma 2.1 を適用する)",
    sectionId := "ncbelyi-lemma-2-2" }

def lemma_2_2.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 3,
    item := "Lemma 2.2",
    sectionId := "ncbelyi-lemma-2-2" }

end ABC3.Found.NCBelyi
