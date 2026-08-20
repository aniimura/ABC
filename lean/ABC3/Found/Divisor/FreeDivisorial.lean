/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import ABC3.Found.FrdI.MonoidVocabulary
import Mathlib.Data.Finsupp.Order

/-!
# 自由可換単系 `ℕ[S]` は divisorial

★★[FrdI] `Example 6.1` は `Φ(L) ⊆ ℤ≥0[D_L]`(台が `D_L` の Cartier 有効因子)を
**divisorial monoid** として使う。その土台が `ℕ[S] = (S →₀ ℕ)` である
(素因子がすべて Cartier なら `Φ(L)` は `ℕ[D_L]` そのもの)。

## ★中身

`IsDivisorial M = (integral ∧ saturated ∧ of characteristic type) ∧ sharp` のうち

* **sharp** —— 成分ごとに `a s + b s = 0` から出る
* **integral** —— `ℕ[S]` は消約的(`isIntegralMonoid_of_isCancelAdd`)
* **of characteristic type** —— sharp から出る(`isOfCharacteristicType_of_isSharp`)
* **saturated** —— ★ここだけ内容がある。`M^gp` の元を差 `x - y` で書き、
  `n·y ≤ n·x`(代数的順序)から `y ≤ x` を出す。`ℕ` では `n > 0` の消約でよい。

★`isSaturatedMonoid_of_cancel_of_nsmul_le` は**一般の消約的単系に使える形**にしてある。
-/

namespace ABC3.Found.FrdI

open AddLocalization

/-- ★★**消約的で「`n·y ≤ n·x` ⟹ `y ≤ x`」(代数的順序)なら saturated**。

★`M^gp` の元を `AddLocalization.ind` で差の形に開き、消約で戻す。 -/
theorem isSaturatedMonoid_of_cancel_of_nsmul_le (M : Type*) [AddCommMonoid M] [IsCancelAdd M]
    (h : ∀ (x y : M) (n : ℕ), 0 < n → (∃ w, n • y + w = n • x) → ∃ w, y + w = x) :
    IsSaturatedMonoid M := by
  intro a n hn ha
  induction a using AddLocalization.ind with
  | _ p =>
    obtain ⟨x, y⟩ := p
    obtain ⟨z, hz⟩ := ha
    rw [toGp, show (n • AddLocalization.mk (x, y).1 (x, y).2)
        = AddLocalization.mk (n • x) (n • y) from mk_nsmul n x y,
      AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists] at hz
    obtain ⟨c, hc⟩ := hz
    have hz' : n • (y : M) + z = n • x := by simpa using hc
    obtain ⟨w, hw⟩ := h x y n hn ⟨z, hz'⟩
    refine ⟨w, ?_⟩
    rw [toGp, AddLocalization.mk_eq_mk_iff, AddLocalization.r_iff_exists]
    exact ⟨0, by simp [← hw, add_comm]⟩

/-- ★★`ℕ[S]` は sharp —— 成分ごとに見るだけ。 -/
theorem isSharp_finsuppNat (S : Type*) : IsSharp (S →₀ ℕ) := by
  intro a ha
  obtain ⟨u, hu⟩ := ha
  have h0 : a + (u.neg : S →₀ ℕ) = 0 := by
    have h := u.val_neg
    rw [hu] at h
    exact h
  ext s
  have h1 := congrArg (fun t : S →₀ ℕ => t s) h0
  simp only [Finsupp.add_apply, Finsupp.coe_zero, Pi.zero_apply] at h1
  simp only [Finsupp.coe_zero, Pi.zero_apply]
  omega

/-- ★★★★**自由可換単系 `ℕ[S]` は divisorial**。

★[FrdI] `Example 6.1` の `Φ(L) ⊆ ℤ≥0[D_L]` が乗る土台である。 -/
theorem isDivisorial_finsuppNat (S : Type*) : IsDivisorial (S →₀ ℕ) := by
  have hsharp := isSharp_finsuppNat S
  refine ⟨⟨isIntegralMonoid_of_isCancelAdd _, ?_, isOfCharacteristicType_of_isSharp _ hsharp⟩,
    hsharp⟩
  refine isSaturatedMonoid_of_cancel_of_nsmul_le _ ?_
  rintro x y n hn ⟨w, hw⟩
  refine ⟨x - y, ?_⟩
  ext s
  have h1 : n * y s + w s = n * x s := by
    have h := congrArg (fun t : S →₀ ℕ => t s) hw
    simpa [Finsupp.add_apply, Finsupp.smul_apply, smul_eq_mul] using h
  have h2 : n * y s ≤ n * x s := by omega
  have h3 : y s ≤ x s := Nat.le_of_mul_le_mul_left h2 hn
  simp only [Finsupp.add_apply, Finsupp.tsub_apply]
  omega

/-! ### ★出典の紐付け -/

/-- ★locator —— `Example 6.1` の `Φ(L) ⊆ ℤ≥0[D_L]` が divisorial であることの土台。 -/
def isDivisorial_finsuppNat.src : ABC3.Meta.Source :=
  { paper := "FrdI", pdfPage := 109,
    item := "Example 6.1 — Φ(L) の土台 ℕ[D_L] は divisorial",
    sectionId := "frdi-example-6-1" }

end ABC3.Found.FrdI
