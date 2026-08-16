import ABC3.Meta.Claim
import Mathlib.Analysis.SpecialFunctions.Pow.Real
import Mathlib.Data.Finset.Max

/-!
# [NCBelyi] Lemma 2.2 —— `λ` による正規化(`Found`)

原典: S. Mochizuki, *Noncritical Belyi Maps* [NCBelyi]、物理 p.4。
**260 dpi 目視確認 2026-08-17**。

原文 (NCBelyi p.4):
> λ. Then, so long as |S| ≥ 4, the polynomial “f(x) + f0” of Lemma 2.1 determines

## ★★原文の「for some appropriate positive rational number λ」を明示する

`Lemma 2.2` の帰納法は、`Lemma 2.1` を **`λ·S`** に適用する。
★原文は `λ` を「適当な正の有理数」としか書かない。**それを明示する。**

`S\{0,∞}` の元を小さい順に `α₁ < α₂ < …` とすると、

> **`λ ≝ 1/α₂`**

と取れば `λ·S` が `Lemma 2.1` の (i)(ii)(iii) を満たす:
- `λα₁ = α₁/α₂ ∈ (0,1)` —— これが `Lemma 2.1` の `r`
- `λα₂ = 1`
- `j ≥ 3` では `λα_j > 1`

★★**`|S| ≥ 4` が要るのは、`α₂` が存在するため**である
(`S` は `0`, `∞` を含むので、`|S| ≥ 4` ⟺ `|S\{0,∞}| ≥ 2`)。
★原文の `|S| ≥ 4` の由来がこれである。

## ★条件 (iv)(`β/α ≥ C`)は `λ` 倍で不変

`(λβ)/(λα) = β/α` なので、正規化しても条件は保たれる。
★これが「`λ·S` に適用してよい」ことの中身である。
-/

namespace ABC3.Found.NCBelyi

/-- ★★**正規化する `λ` の存在**。

正の有理数の有限集合 `T`(= `S\{0,∞}`)が 2 元以上を持てば、
`λ > 0` と `r ∈ (0,1)` があって `λ·T` の元は `r`、`1`、または `> 1` のいずれかである。

★★`λ ≝ 1/α₂`(`α₂` は 2 番目に小さい元)である。 -/
theorem exists_normalizing_scale (T : Finset ℚ) (hpos : ∀ α ∈ T, 0 < α)
    (h2 : 2 ≤ T.card) :
    ∃ lam r : ℚ, 0 < lam ∧ 0 < r ∧ r < 1 ∧
      (1 : ℚ) ∈ T.image (fun t => lam * t) ∧
      r ∈ T.image (fun t => lam * t) ∧
      ∀ α ∈ T.image (fun t => lam * t), α = r ∨ α = 1 ∨ 1 < α := by
  classical
  have hne : T.Nonempty := Finset.card_pos.1 (by omega)
  set a1 : ℚ := T.min' hne with ha1
  have ha1T : a1 ∈ T := T.min'_mem hne
  have hT' : (T.erase a1).Nonempty := by
    rw [← Finset.card_pos, Finset.card_erase_of_mem ha1T]
    omega
  set a2 : ℚ := (T.erase a1).min' hT' with ha2
  have ha2T' : a2 ∈ T.erase a1 := (T.erase a1).min'_mem hT'
  have ha2T : a2 ∈ T := Finset.mem_of_mem_erase ha2T'
  have ha2ne : a2 ≠ a1 := Finset.ne_of_mem_erase ha2T'
  have ha1a2 : a1 < a2 := lt_of_le_of_ne (T.min'_le a2 ha2T) (Ne.symm ha2ne)
  have ha20 : 0 < a2 := hpos a2 ha2T
  have ha10 : 0 < a1 := hpos a1 ha1T
  refine ⟨1 / a2, a1 / a2, div_pos one_pos ha20, div_pos ha10 ha20, ?_, ?_, ?_, ?_⟩
  · rw [div_lt_one ha20]; exact ha1a2
  · exact Finset.mem_image.2 ⟨a2, ha2T, by field_simp⟩
  · exact Finset.mem_image.2 ⟨a1, ha1T, by ring⟩
  · intro α hα
    obtain ⟨t, htT, rfl⟩ := Finset.mem_image.1 hα
    by_cases h1 : t = a1
    · left; rw [h1]; ring
    by_cases h2' : t = a2
    · right; left; rw [h2']; field_simp
    · right; right
      have htT' : t ∈ T.erase a1 := Finset.mem_erase.2 ⟨h1, htT⟩
      have hle : a2 ≤ t := (T.erase a1).min'_le t htT'
      have hlt : a2 < t := lt_of_le_of_ne hle (Ne.symm h2')
      rw [one_div, inv_mul_eq_div, lt_div_iff₀ ha20, one_mul]
      exact hlt

/-- ★**条件 `β/α ≥ C` は `λ` 倍で不変**。

★これが「`λ·S` に `Lemma 2.1` を適用してよい」ことの中身である。 -/
theorem ratio_scale_invariant {lam β α : ℚ} (hlam : 0 < lam) (hα : 0 < α) {C : ℚ}
    (h : C ≤ β / α) : C ≤ (lam * β) / (lam * α) := by
  rw [mul_div_mul_left _ _ hlam.ne']
  exact h

/-! ## ★`|S′| < |S|` —— 2 点が同じ値へ行けば像は小さくなる -/

/-- ★**衝突があれば像の要素数は減る**。

原文 (NCBelyi p.4):
> which the cardinalities of S , S satisfy |S | < |S|.

★`Lemma 2.1` の `f` は `f(0) = f(1) = 0` なので、
`0` と `1` が**同じ値へ行く**。ここから `|S′| < |S|` が出る。
★★原文は `|S′| < |S|` とだけ書くが、**その根拠は `f(0) = f(1)` である**。 -/
theorem card_image_lt_of_collision {α β : Type*} [DecidableEq α] [DecidableEq β]
    {T : Finset α} {g : α → β} {x y : α}
    (hx : x ∈ T) (hy : y ∈ T) (hxy : x ≠ y) (hg : g x = g y) :
    (T.image g).card < T.card := by
  classical
  have hsub : T.image g = (T.erase x).image g := by
    refine Finset.Subset.antisymm (fun z hz => ?_) ?_
    · obtain ⟨t, htT, rfl⟩ := Finset.mem_image.1 hz
      by_cases h : t = x
      · subst h
        exact Finset.mem_image.2 ⟨y, Finset.mem_erase.2 ⟨Ne.symm hxy, hy⟩, hg.symm⟩
      · exact Finset.mem_image.2 ⟨t, Finset.mem_erase.2 ⟨h, htT⟩, rfl⟩
    · exact Finset.image_subset_image (Finset.erase_subset _ _)
  rw [hsub]
  calc ((T.erase x).image g).card ≤ (T.erase x).card := Finset.card_image_le
    _ < T.card := Finset.card_erase_lt_of_mem hx

/-! ## ★出典の紐付け(`.src`) -/

def exists_normalizing_scale.src : ABC3.Meta.Source :=
  { paper := "NCBelyi", pdfPage := 3,
    item := "Lemma 2.2(λ による正規化のみ)",
    sectionId := "ncbelyi-lemma-2-2" }

end ABC3.Found.NCBelyi
