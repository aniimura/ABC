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

end ABC3.Found.NCBelyi
