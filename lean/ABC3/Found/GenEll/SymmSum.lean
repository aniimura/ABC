/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Algebra.BigOperators.Intervals
import ABC3.Meta.Claim

/-!
# 第 956 ブロック —— **★★★★★★★★★★★★`i ↦ l-i` で対称な和は 2 で割れる**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——(D3) の (d) の鍵

Vélu の `w` は `2 * w = ∑ (…)` という形で受けている（`2` で割らないため）。
★だから `w` を作るには**右辺の和が 2 で割れる**ことが要る。

☆`μ_l` の点は `i ↦ l-i`（すなわち反転）で対になる。
`l` は奇数なのでこの対応に**不動点がない**——`i = l-i` は `2i = l` を意味するから。
★したがって添字集合はきちんと 2 つずつの組に分かれ、
対称な関数の和は `2 • (半分の和)` になる。

☆本ブロックはこれを**楕円曲線とは無関係に**、純粋な `Finset` の事実として取る。

| 定理 | 内容 |
|---|---|
| `erase_zero_range_eq_union` | ★`(range l).erase 0` は前半と後半の素分割 |
| `sum_eq_two_nsmul_of_symm` | ★★★★★★★★★★★★**対称な和は `2 • (半分)`** |
| `exists_two_mul_of_symm` | ★★★★★★★★**環での `2 * w` の形** |
-/

namespace ABC3.Found.GenEll

open Finset

section SymmSum

variable {l : ℕ}

/-- ★**`(range l).erase 0` は前半 `Icc 1 m` とその反転像の素分割**。

☆`l = 2m+1`（奇数）のとき、`1 ≤ i ≤ l-1` は
`i ≤ m` か `i ≥ m+1` のどちらかであり、後者は `l-i ≤ m` に対応する。 -/
theorem erase_zero_range_eq_union (m : ℕ) :
    (range (2 * m + 1)).erase 0
      = (Icc 1 m) ∪ ((Icc 1 m).image (fun i => 2 * m + 1 - i)) := by
  ext i
  simp only [mem_erase, mem_range, mem_union, mem_Icc, mem_image]
  constructor
  · rintro ⟨h0, hi⟩
    by_cases hm : i ≤ m
    · exact Or.inl ⟨by omega, hm⟩
    · exact Or.inr ⟨2 * m + 1 - i, ⟨by omega, by omega⟩, by omega⟩
  · rintro (⟨h1, h2⟩ | ⟨j, ⟨hj1, hj2⟩, rfl⟩)
    · exact ⟨by omega, by omega⟩
    · exact ⟨by omega, by omega⟩

/-- ★前半とその反転像は交わらない。 -/
theorem disjoint_Icc_image_sub (m : ℕ) :
    Disjoint (Icc 1 m) ((Icc 1 m).image (fun i => 2 * m + 1 - i)) := by
  rw [Finset.disjoint_left]
  intro a ha hb
  simp only [mem_Icc] at ha
  simp only [mem_image, mem_Icc] at hb
  obtain ⟨j, ⟨hj1, hj2⟩, rfl⟩ := hb
  omega

/-- ★★★★★★★★★★★★**`i ↦ l-i` で対称な和は `2 • (半分の和)`**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★`l` が奇数であることが本質である——不動点がないからちょうど 2 つずつの組になる。 -/
theorem sum_eq_two_nsmul_of_symm {M : Type} [AddCommMonoid M] (m : ℕ) (f : ℕ → M)
    (hf : ∀ i ∈ Icc 1 m, f (2 * m + 1 - i) = f i) :
    ∑ i ∈ (range (2 * m + 1)).erase 0, f i = 2 • (∑ i ∈ Icc 1 m, f i) := by
  rw [erase_zero_range_eq_union m, Finset.sum_union (disjoint_Icc_image_sub m)]
  have hinj : ∀ a ∈ (Icc 1 m : Finset ℕ), ∀ b ∈ (Icc 1 m : Finset ℕ),
      2 * m + 1 - a = 2 * m + 1 - b → a = b := by
    intro a ha b hb hab
    simp only [mem_Icc] at ha hb
    omega
  rw [Finset.sum_image hinj]
  rw [Finset.sum_congr rfl hf]
  rw [two_nsmul]

/-- ★★★★★★★★**環での `2 * w` の形**——Vélu の `w` を作る形である。

★★★★**2026-09-01（第 956）**——`hw : 2 * w = ∑ (…)` を満たす `w` は
**和が `i ↦ l-i` で対称でありさえすれば作れる**。 -/
theorem exists_two_mul_of_symm {A : Type} [CommRing A] (m : ℕ) (f : ℕ → A)
    (hf : ∀ i ∈ Icc 1 m, f (2 * m + 1 - i) = f i) :
    ∃ w : A, 2 * w = ∑ i ∈ (range (2 * m + 1)).erase 0, f i := by
  refine ⟨∑ i ∈ Icc 1 m, f i, ?_⟩
  rw [sum_eq_two_nsmul_of_symm m f hf, two_nsmul]
  exact two_mul _

end SymmSum

/-! ## ★出典の紐付け(`.src`) -/

def sum_eq_two_nsmul_of_symm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(i ↦ l-i で対称な和は 2 • (半分)。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_two_mul_of_symm.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(Vélu の w は和が対称なら作れる。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
