import ABC3.Found.NCBelyi.Normalize
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

/-! ## ★出典の紐付け(`.src`) -/

def lemma_2_2_base.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 4,
    item := "Lemma 2.2(基底段——原文が書いていない |S| ≤ 3 の場合)",
    sectionId := "ncbelyi-lemma-2-2" }

end ABC3.Found.NCBelyi
