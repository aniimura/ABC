/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.Algebra.BigOperators.Intervals
import Mathlib.Algebra.BigOperators.Ring.Finset
import ABC3.Meta.Claim
import ABC3.Found.GenEll.Velu

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

/-- ★★★★★★★★★★★★★★★★**対ごとに偶なら和も偶**——
`sum_eq_two_nsmul_of_symm` の一般形。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 957）の訂正**——Vélu の和の項は
**反転で不変ではない**。Tate 曲線は `a₁ = 1` なので

    `veluV2 = 3x² + a₄ - y`

であり、`y ↦ -y-x` で変わる。☆変わらないのは `veluU = (2y+x)²` の方だけである。

★だが**対の和は偶である**:

    `g(i) + g(l-i) = 2·veluU_i + 2·x_i·(veluV2_i + veluV2_{l-i})`

なので、対称性ではなく**対ごとの偶性**を受けるのが正しい。
☆証人 `c` を関数として受けるので選択公理も要らない。 -/
theorem two_mul_sum_eq_of_pair_even {A : Type} [CommRing A] (m : ℕ) (g c : ℕ → A)
    (h : ∀ i ∈ Icc 1 m, g i + g (2 * m + 1 - i) = 2 * c i) :
    2 * (∑ i ∈ Icc 1 m, c i) = ∑ i ∈ (range (2 * m + 1)).erase 0, g i := by
  have hinj : ∀ a ∈ (Icc 1 m : Finset ℕ), ∀ b ∈ (Icc 1 m : Finset ℕ),
      2 * m + 1 - a = 2 * m + 1 - b → a = b := by
    intro a ha b hb hab
    simp only [mem_Icc] at ha hb
    omega
  rw [erase_zero_range_eq_union m, Finset.sum_union (disjoint_Icc_image_sub m),
    Finset.sum_image hinj, ← Finset.sum_add_distrib, Finset.sum_congr rfl h,
    Finset.mul_sum]

/-- ★★★★★★★★★★★★**対ごとに偶なら `w` が作れる**。

★これが Vélu の `hw : 2 * w = ∑ (…)` を満たす `w` の作り方である。 -/
theorem exists_two_mul_of_pair_even {A : Type} [CommRing A] (m : ℕ) (g c : ℕ → A)
    (h : ∀ i ∈ Icc 1 m, g i + g (2 * m + 1 - i) = 2 * c i) :
    ∃ w : A, 2 * w = ∑ i ∈ (range (2 * m + 1)).erase 0, g i :=
  ⟨∑ i ∈ Icc 1 m, c i, two_mul_sum_eq_of_pair_even m g c h⟩

/-! ## ★★★★★★★★★★★★★★★★(d3)——Vélu の `w` を作る -/

open WeierstrassCurve in
/-- ★★★★★★★★★★★★★★★★**添字の反転が点の反転なら、Vélu の `w` は環の中で作れる**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 960）**——これが (D3) の (d3) である。
☆三つを並べるだけである:

1. `two_mul_sum_eq_of_pair_even`（第 957）——対ごとに偶なら和も偶
2. `veluTerm_pair_even`（第 958）——Vélu の項の対の和は偶
3. 仮説 `hX`・`hY`——添字の反転が点の反転（`tateXpair_mu_inv`、第 959）

★`l = 2m+1`（奇数）であることが本質である——不動点がないから対になる。 -/
theorem exists_veluW_of_inv {A : Type} [CommRing A] (W : WeierstrassCurve A) (m : ℕ)
    (X Y : ℕ → A)
    (hX : ∀ i ∈ Icc 1 m, X (2 * m + 1 - i) = X i)
    (hY : ∀ i ∈ Icc 1 m, Y (2 * m + 1 - i) = W.toAffine.negY (X i) (Y i)) :
    ∃ w : A, 2 * w = ∑ i ∈ (range (2 * m + 1)).erase 0,
      (veluU W (X i) (Y i) + 2 * (veluV2 W (X i) (Y i) * X i)) := by
  refine exists_two_mul_of_pair_even m
    (fun i => veluU W (X i) (Y i) + 2 * (veluV2 W (X i) (Y i) * X i))
    (fun i => veluU W (X i) (Y i)
      + X i * (veluV2 W (X i) (Y i)
          + veluV2 W (X i) (W.toAffine.negY (X i) (Y i)))) ?_
  intro i hi
  rw [hX i hi, hY i hi]
  exact veluTerm_pair_even W (X i) (Y i)

/-! ## ★★★★★★★★第 1149 —— `l = 2` のときの `w` -/

open WeierstrassCurve in
/-- ★★★★★★★★★★★★**`l = 2` でも Vélu の `w` は環の中で作れる**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 1149）**——`exists_veluW_of_inv`（第 960）は
`l = 2m+1`（奇）を使っている。☆`l = 2` では添字集合が
`(range 2).erase 0 = {1}` と一点だけになり、対にならない。

★だが**その一点は 2-捻れ**である——`Y 1 = negY (X 1) (Y 1)`。
すなわち `veluGy = -2y - a₁x - a₃ = 0` なので `veluU = veluGy² = 0` となり、

    `∑ (veluU + 2·veluV2·x) = 2·(veluV2·x)`

である。☆つまり `w = veluV2·x` がそのまま取れる——**割り算をしない**。 -/
theorem exists_veluW_two {A : Type} [CommRing A] (W : WeierstrassCurve A) (X Y : ℕ → A)
    (hY : Y 1 = W.toAffine.negY (X 1) (Y 1)) :
    ∃ w : A, 2 * w = ∑ i ∈ (range 2).erase 0,
      (veluU W (X i) (Y i) + 2 * (veluV2 W (X i) (Y i) * X i)) := by
  refine ⟨veluV2 W (X 1) (Y 1) * X 1, ?_⟩
  have hs : (range 2).erase 0 = {1} := by decide
  rw [hs, Finset.sum_singleton]
  have hgy : veluGy W (X 1) (Y 1) = 0 := by
    rw [WeierstrassCurve.Affine.negY] at hY
    simp only [veluGy]
    linear_combination -hY
  rw [veluU, hgy]
  ring

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


def two_mul_sum_eq_of_pair_even.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(対ごとに偶なら和も偶。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_two_mul_of_pair_even.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(対ごとに偶なら Vélu の w が作れる。★無条件)",
    sectionId := "genell-lemma-3-5" }


def exists_veluW_of_inv.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(添字の反転が点の反転なら Vélu の w は環の中で作れる。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_veluW_two.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(l = 2 でも Vélu の w は作れる——一点は 2-捻れで veluU = 0。★無条件)",
    sectionId := "genell-lemma-3-5" }

def exists_veluW_two.needs : List ABC3.Meta.ProofObligation :=
  [ .implicitStep
      ("★★★★**2026-09-01（第 1149）の測定**——`exists_veluW_of_inv` は " ++
       "`l = 2m+1`（奇）を使うが、`l = 2` では添字集合が `{1}` と一点だけで対にならない。" ++
       "☆しかしその一点は 2-捻れなので `veluGy = 0`、したがって `veluU = 0` であり、" ++
       "`∑ (veluU + 2·veluV2·x) = 2·(veluV2·x)` となる。★割り算をしない。") 6 ]

end ABC3.Found.GenEll
