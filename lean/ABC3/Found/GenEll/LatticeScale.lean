import ABC3.Found.GenEll.LatticeHalfPeriod

/-!
# GenEll 第 333 ブロック —— **★★★★★★束のスケール変換と判別式**(`Found`)

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any

## ★★★★★★★これで判別式の非消失が「正規化された束」に帰着する

G8 の `htFalt` を固定する連鎖の (i)「束の判別式 `g₂³ − 27g₃² ≠ 0`」に向けた 2 段目。

★束を `c` 倍すると Eisenstein 級数は `c⁻ⁿ` 倍される:

    G(cΛ, n) = c⁻ⁿ · G(Λ, n)

★★したがって `g₂ = 60G₄` は `c⁻⁴` 倍、`g₃ = 140G₆` は `c⁻⁶` 倍、
**判別式 `g₂³ − 27g₃²` は `c⁻¹²` 倍**される。
★★★**非消失はスケールに依らない**(`latticeDisc_ne_zero_iff`)ので、
`Λ = ω₁·(ℤ + τℤ)` と書けば `Λ_τ = ℤ + τℤ` の場合に帰着する。

## ★★★★★★これが効く先——モジュラー判別式へ

正規化された束では、古典的に

    g₂(Λ_τ)³ − 27g₃(Λ_τ)² = (2π)¹² · Δ(τ),    Δ = η²⁴

であり、★**mathlib は `Δ = η²⁴` とその非消失を持つ**(`NumberTheory/ModularForms/Discriminant.lean`)。
★★したがって (i) は「格子 Eisenstein 級数 `G` とモジュラー Eisenstein 級数 `E` の
正規化の対応」に帰着する——**偏角の原理は要らない**。

★★★★★これは当初見込んでいた道(基本平行四辺形の周での偏角の原理)より**軽い**。

## ★★本ブロックで取れるもの

| 定理 | 内容 |
|---|---|
| `scalePair` | ★束の `c` 倍 |
| `scalePair_lattice` | ★★束は `c` 倍で写る |
| `G_scalePair` | ★★★★★**`G(cΛ, n) = c⁻ⁿ G(Λ, n)`** |
| `latticeDisc` | ★判別式 `g₂³ − 27g₃²` |
| `g₂_scalePair`・`g₃_scalePair` | ★★★`c⁻⁴`・`c⁻⁶` 倍 |
| `latticeDisc_scalePair` | ★★★★★★**判別式は `c⁻¹²` 倍** |
| `latticeDisc_ne_zero_iff` | ★★★★★**非消失はスケールに依らない** |
-/

namespace ABC3.Found.GenEll

open PeriodPair

/-! ## ★束の `c` 倍 -/

/-- ★**束の `c` 倍**(`c ≠ 0`)。 -/
noncomputable def scalePair (L : PeriodPair) (c : ℂ) (hc : c ≠ 0) : PeriodPair where
  ω₁ := c * L.ω₁
  ω₂ := c * L.ω₂
  indep := by
    have hker : LinearMap.ker (LinearMap.mulLeft ℝ c) = ⊥ := by
      rw [LinearMap.ker_eq_bot]
      intro x y hxy
      simpa [LinearMap.mulLeft_apply, mul_right_inj' hc] using hxy
    have h := L.indep.map' (LinearMap.mulLeft ℝ c) hker
    have heq : (⇑(LinearMap.mulLeft ℝ c) ∘ ![L.ω₁, L.ω₂]) = ![c * L.ω₁, c * L.ω₂] := by
      funext i
      fin_cases i <;> simp
    rwa [heq] at h

/-- ★★束も `c` 倍で写る。 -/
theorem scalePair_lattice (L : PeriodPair) (c : ℂ) (hc : c ≠ 0) :
    (scalePair L c hc).lattice = L.lattice.map (LinearMap.mulLeft ℤ c) := by
  rw [PeriodPair.lattice, PeriodPair.lattice, Submodule.map_span]
  congr 1
  ext x
  simp only [Set.mem_image, Set.mem_insert_iff, Set.mem_singleton_iff, LinearMap.mulLeft_apply]
  constructor
  · rintro (rfl | rfl)
    · exact ⟨L.ω₁, Or.inl rfl, rfl⟩
    · exact ⟨L.ω₂, Or.inr rfl, rfl⟩
  · rintro ⟨y, (rfl | rfl), rfl⟩
    · exact Or.inl rfl
    · exact Or.inr rfl

/-! ## ★★★★★Eisenstein 級数のスケール則 -/

set_option maxHeartbeats 1600000 in
/-- ★★★★★**`G(cΛ, n) = c⁻ⁿ · G(Λ, n)`**——和の添字を `l ↦ cl` で取り替える。 -/
theorem G_scalePair (L : PeriodPair) (c : ℂ) (hc : c ≠ 0) (n : ℕ) :
    (scalePair L c hc).G n = (c ^ n)⁻¹ * L.G n := by
  have hlat := scalePair_lattice L c hc
  have hmem : ∀ l : L.lattice, (c * (l : ℂ)) ∈ (scalePair L c hc).lattice := by
    intro l
    rw [hlat]
    exact ⟨l, l.2, rfl⟩
  have hmem' : ∀ l' : (scalePair L c hc).lattice, ((l' : ℂ) / c) ∈ L.lattice := by
    intro l'
    have h2 : (l' : ℂ) ∈ L.lattice.map (LinearMap.mulLeft ℤ c) := by
      rw [← hlat]; exact l'.2
    obtain ⟨y, hy, hyeq⟩ := h2
    have h3 : (l' : ℂ) = c * y := hyeq.symm
    have h4 : c * y / c = y := by field_simp
    rw [h3, h4]
    exact hy
  let e : L.lattice ≃ (scalePair L c hc).lattice :=
    { toFun := fun l => ⟨c * l, hmem l⟩
      invFun := fun l' => ⟨(l' : ℂ) / c, hmem' l'⟩
      left_inv := by
        intro l
        ext
        show c * (l : ℂ) / c = (l : ℂ)
        field_simp
      right_inv := by
        intro l'
        ext
        show c * ((l' : ℂ) / c) = (l' : ℂ)
        field_simp }
  rw [PeriodPair.G, PeriodPair.G,
    ← e.tsum_eq (fun l' : (scalePair L c hc).lattice => (((l' : ℂ)) ^ n)⁻¹), ← tsum_mul_left]
  congr 1
  funext l
  show ((c * (l : ℂ)) ^ n)⁻¹ = (c ^ n)⁻¹ * ((l : ℂ) ^ n)⁻¹
  rw [mul_pow, mul_inv]

/-! ## ★★★★★★判別式のスケール則 -/

/-- ★★★★★束の判別式 `g₂³ − 27g₃²`。 -/
noncomputable def latticeDisc (L : PeriodPair) : ℂ := L.g₂ ^ 3 - 27 * L.g₃ ^ 2

/-- ★★★`g₂` は `c⁻⁴` 倍。 -/
theorem g₂_scalePair (L : PeriodPair) (c : ℂ) (hc : c ≠ 0) :
    (scalePair L c hc).g₂ = (c ^ 4)⁻¹ * L.g₂ := by
  rw [PeriodPair.g₂, PeriodPair.g₂, G_scalePair]; ring

/-- ★★★`g₃` は `c⁻⁶` 倍。 -/
theorem g₃_scalePair (L : PeriodPair) (c : ℂ) (hc : c ≠ 0) :
    (scalePair L c hc).g₃ = (c ^ 6)⁻¹ * L.g₃ := by
  rw [PeriodPair.g₃, PeriodPair.g₃, G_scalePair]; ring

/-- ★★★★★★**判別式は `c⁻¹²` 倍される**。

原文 (GenEll p.17):
> Proposition 3.4. (Faltings Heights and the Divisor at Infinity) For any -/
theorem latticeDisc_scalePair (L : PeriodPair) (c : ℂ) (hc : c ≠ 0) :
    latticeDisc (scalePair L c hc) = (c ^ 12)⁻¹ * latticeDisc L := by
  rw [latticeDisc, latticeDisc, g₂_scalePair, g₃_scalePair]
  field_simp

/-- ★★★★★**判別式の非消失はスケールに依らない**——正規化された束に帰着できる。 -/
theorem latticeDisc_ne_zero_iff (L : PeriodPair) (c : ℂ) (hc : c ≠ 0) :
    latticeDisc (scalePair L c hc) ≠ 0 ↔ latticeDisc L ≠ 0 := by
  rw [latticeDisc_scalePair]
  constructor
  · intro h h0
    exact h (by rw [h0, mul_zero])
  · intro h h0
    exact h ((mul_eq_zero.1 h0).resolve_left (inv_ne_zero (pow_ne_zero 12 hc)))

/-! ## ★出典の紐付け(`.src`) -/

def G_scalePair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

def latticeDisc_scalePair.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Proposition 3.4(Faltings Heights and the Divisor at Infinity)",
    sectionId := "genell-prop-3-4" }

end ABC3.Found.GenEll
