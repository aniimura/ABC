/-
Copyright (c) 2026 ABC3. All rights reserved.
-/
import Mathlib.RingTheory.RootsOfUnity.Lemmas
import ABC3.Meta.Claim

/-!
# 第 951 ブロック —— **★★★★★★★★★★★★`l` が単元なら `1 - ζ^i` も単元**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★これは何か——仮説を 1 つ減らす

`tateParam_quot_velu_dvr`（第 927）は `hlu : IsUnit (l : R)` と
`hu : ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ^i)` の**両方**を受けていた。

★だが分円分解の古典的な等式

    `∏_{i=1}^{l-1} (1 - ζ^i) = l`

があるので、**`hu` は `hlu` から出る**。
☆積が単元なら各因子も単元だからである（可換環）。

★mathlib の `IsPrimitiveRoot.prod_one_sub_pow_eq_order` は
`∏_{k ∈ range n} (1 - μ^{k+1}) = n+1` の形であるので、
本ブロックはそれを `(range l).erase 0` の形に乗せ換える。

| 定理 | 内容 |
|---|---|
| `range_erase_zero_eq_image_succ` | ★`(range l).erase 0 = (range (l-1)).image (·+1)` |
| `prod_one_sub_pow_erase` | ★★★★`∏_{i=1}^{l-1} (1 - ζ^i) = l` |
| `isUnit_one_sub_pow_of_isUnit_natCast` | ★★★★★★★★★★★★**`hu` は `hlu` から出る** |
-/

namespace ABC3.Found.GenEll

open Finset

/-- ★**`(range (n+1)).erase 0` は `(range n)` の `(·+1)` による像**。 -/
theorem range_erase_zero_eq_image_succ (n : ℕ) :
    (range (n + 1)).erase 0 = (range n).image (fun k => k + 1) := by
  ext i
  simp only [mem_erase, mem_range, mem_image]
  constructor
  · rintro ⟨h0, hi⟩
    exact ⟨i - 1, by omega, by omega⟩
  · rintro ⟨k, hk, rfl⟩
    omega

/-- ★★★★**`∏_{i=1}^{l-1} (1 - ζ^i) = l`**——分円分解の古典的な等式。

☆mathlib の `IsPrimitiveRoot.prod_one_sub_pow_eq_order` を
`(range l).erase 0` の形に乗せ換えただけである。 -/
theorem prod_one_sub_pow_erase {R : Type} [CommRing R] [IsDomain R] {l : ℕ} (hl : 0 < l)
    {ζ : R} (hζ : IsPrimitiveRoot ζ l) :
    ∏ i ∈ (range l).erase 0, (1 - ζ ^ i) = (l : R) := by
  obtain ⟨n, rfl⟩ : ∃ n, l = n + 1 := ⟨l - 1, by omega⟩
  rw [range_erase_zero_eq_image_succ n,
    Finset.prod_image (fun a _ b _ h => by omega)]
  rw [hζ.prod_one_sub_pow_eq_order]
  push_cast
  ring

/-- ★★★★★★★★★★★★**`l` が単元なら `1 - ζ^i` も単元**。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★★**2026-09-01（第 951）**——これで `tateParam_quot_velu_dvr`（第 927）の
**仮説 `hu` が 1 つ消える**。積が `l` に等しいので、
積が単元なら各因子も単元である。 -/
theorem isUnit_one_sub_pow_of_isUnit_natCast {R : Type} [CommRing R] [IsDomain R]
    {l : ℕ} (hl : 0 < l) {ζ : R} (hζ : IsPrimitiveRoot ζ l) (hlu : IsUnit ((l : R))) :
    ∀ i ∈ (range l).erase 0, IsUnit (1 - ζ ^ i) := by
  intro i hi
  rw [← prod_one_sub_pow_erase hl hζ, ← Finset.mul_prod_erase _ _ hi] at hlu
  exact isUnit_of_mul_isUnit_left hlu

/-! ## ★出典の紐付け(`.src`) -/

def prod_one_sub_pow_erase.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(∏_{i=1}^{l-1} (1 - ζ^i) = l。★無条件)",
    sectionId := "genell-lemma-3-5" }

def isUnit_one_sub_pow_of_isUnit_natCast.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(l が単元なら 1 - ζ^i も単元。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
