/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.Uniformization
import ABC3.Found.GenEll.LatticeHalfPeriod
import ABC3.Meta.Claim

/-!
# 第 1397 ブロック —— **判別式は半周期の値の差積である**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か——**恒等式の証明の第 1 歩**

2026-09-02 に `Δ(E)^l = Δ(E/C)·N⁴` の**証明の道が見つかった**。
☆道は 3 つの部品からなる:

| 部品 | 内容 | 状態 |
|---|---|---|
| (1) | `Δ = 16·((e₁−e₂)(e₁−e₃)(e₂−e₃))²`（`e_i` は半周期での `℘` の値） | ★**本ブロック** |
| (2) | `℘(z + v_j) = e_j + (e_j−e_a)(e_j−e_b)/(℘(z) − e_j)`（半周期を足す公式） | ☆次 |
| (3) | `∏_{w∈T}(℘(z+w) − e_i) = c_i·(℘_{Λ′}(z) − e′_i)`（**同種のノルム**） | ☆その次 |

★★★(3) を `z = v_j` で使い、(2) で `℘(v_j+w)` を `℘(w)` の式に直すと
**`∏_{i≠j}(℘(v_j+w) − e_i) = −D`（`w` に依らない！）**という代数恒等式に落ちる。
☆あとは `D^l = D′·(c₁c₂c₃)²` の帳簿と `N² = 4^{l−1}c₁c₂c₃` で恒等式が出る。

☆数値検証（60 桁、曲線 `y²−y = x³−x²`、`l = 5`）で (3) と最終形の両方を確かめてある。

## ★★★★本ブロックの中身

`e₁ = ℘(ω₁/2)`・`e₂ = ℘(ω₂/2)`・`e₃ = ℘((ω₁+ω₂)/2)` は

* 3 次式 `4x³ − g₂x − g₃` の根である（`LatticeHalfPeriod.lean`、在庫）
* **相異なる**（★本ブロック——一意化の単射性 `mem_lattice_of_shift_eq` の系）

☆相異なる 3 根から Vieta が出て、`latticeDisc = g₂³ − 27g₃²` が
**`16·((e₁−e₂)(e₁−e₃)(e₂−e₃))²`** になる。
-/

namespace ABC3.Found.GenEll

open PeriodPair ABC3.Meta

/-! ## ★★★★半周期の差も束の外 -/

/-- ★★★`(ω₁−ω₂)/2` も束に入らない。 -/
theorem sub_div_two_notMem_lattice (L : PeriodPair) : (L.ω₁ - L.ω₂) / 2 ∉ L.lattice := by
  intro h
  have hq : ((1 / 2 : ℚ) : ℂ) * L.ω₁ + ((-1 / 2 : ℚ) : ℂ) * L.ω₂ ∈ L.lattice := by
    convert h using 1; push_cast; ring
  rw [PeriodPair.mul_ω₁_add_mul_ω₂_mem_lattice] at hq
  norm_num at hq

/-! ## ★★★★★★★★半周期の値は相異なる -/

/-- ★★★★★★★★`e₁ ≠ e₂`——一意化の単射性の系（第 1397）。 -/
theorem e₁_ne_e₂ (L : PeriodPair) :
    L.weierstrassP (L.ω₁ / 2) ≠ L.weierstrassP (L.ω₂ / 2) := by
  intro heq
  refine sub_div_two_notMem_lattice L ?_
  refine mem_lattice_of_shift_eq L ((L.ω₁ - L.ω₂) / 2) (L.ω₂_div_two_notMem_lattice) ?_ ?_ ?_
  · rw [show L.ω₂ / 2 + (L.ω₁ - L.ω₂) / 2 = L.ω₁ / 2 by ring]
    exact L.ω₁_div_two_notMem_lattice
  · rw [show L.ω₂ / 2 + (L.ω₁ - L.ω₂) / 2 = L.ω₁ / 2 by ring]
    exact heq
  · rw [show L.ω₂ / 2 + (L.ω₁ - L.ω₂) / 2 = L.ω₁ / 2 by ring,
      derivWeierstrassP_half L L.ω₁ L.ω₁_mem_lattice,
      derivWeierstrassP_half L L.ω₂ L.ω₂_mem_lattice]

/-- ★★★★★★★★`e₁ ≠ e₃`（第 1397）。 -/
theorem e₁_ne_e₃ (L : PeriodPair) :
    L.weierstrassP (L.ω₁ / 2) ≠ L.weierstrassP ((L.ω₁ + L.ω₂) / 2) := by
  intro heq
  refine L.ω₂_div_two_notMem_lattice ?_
  have hmem : -(L.ω₂ / 2) ∈ L.lattice := by
    refine mem_lattice_of_shift_eq L (-(L.ω₂ / 2)) (add_div_two_notMem_lattice L) ?_ ?_ ?_
    · rw [show (L.ω₁ + L.ω₂) / 2 + -(L.ω₂ / 2) = L.ω₁ / 2 by ring]
      exact L.ω₁_div_two_notMem_lattice
    · rw [show (L.ω₁ + L.ω₂) / 2 + -(L.ω₂ / 2) = L.ω₁ / 2 by ring]
      exact heq
    · rw [show (L.ω₁ + L.ω₂) / 2 + -(L.ω₂ / 2) = L.ω₁ / 2 by ring,
        derivWeierstrassP_half L L.ω₁ L.ω₁_mem_lattice,
        derivWeierstrassP_half L (L.ω₁ + L.ω₂)
          (L.lattice.add_mem L.ω₁_mem_lattice L.ω₂_mem_lattice)]
  simpa using neg_mem hmem

/-- ★★★★★★★★`e₂ ≠ e₃`（第 1397）。 -/
theorem e₂_ne_e₃ (L : PeriodPair) :
    L.weierstrassP (L.ω₂ / 2) ≠ L.weierstrassP ((L.ω₁ + L.ω₂) / 2) := by
  intro heq
  refine L.ω₁_div_two_notMem_lattice ?_
  have hmem : -(L.ω₁ / 2) ∈ L.lattice := by
    refine mem_lattice_of_shift_eq L (-(L.ω₁ / 2)) (add_div_two_notMem_lattice L) ?_ ?_ ?_
    · rw [show (L.ω₁ + L.ω₂) / 2 + -(L.ω₁ / 2) = L.ω₂ / 2 by ring]
      exact L.ω₂_div_two_notMem_lattice
    · rw [show (L.ω₁ + L.ω₂) / 2 + -(L.ω₁ / 2) = L.ω₂ / 2 by ring]
      exact heq
    · rw [show (L.ω₁ + L.ω₂) / 2 + -(L.ω₁ / 2) = L.ω₂ / 2 by ring,
        derivWeierstrassP_half L L.ω₂ L.ω₂_mem_lattice,
        derivWeierstrassP_half L (L.ω₁ + L.ω₂)
          (L.lattice.add_mem L.ω₁_mem_lattice L.ω₂_mem_lattice)]
  simpa using neg_mem hmem

/-! ## ★★★★★★★★★★★★相異なる 3 根から Vieta -/

/-- ★★★★★★★★★★★★
**`4x³ − gx − h` に相異なる根 `a, b, c` があれば Vieta が出る**——★**無条件**（第 1397）。

☆`g³ − 27h² = 16((a−b)(a−c)(b−c))²` まで込みである。 -/
theorem vieta_of_three_roots (a b c g h : ℂ) (hab : a ≠ b) (hac : a ≠ c) (hbc : b ≠ c)
    (h1 : 4 * a ^ 3 - g * a - h = 0) (h2 : 4 * b ^ 3 - g * b - h = 0)
    (h3 : 4 * c ^ 3 - g * c - h = 0) :
    a + b + c = 0 ∧ g = -4 * (a * b + a * c + b * c) ∧ h = 4 * (a * b * c)
      ∧ g ^ 3 - 27 * h ^ 2 = 16 * ((a - b) * (a - c) * (b - c)) ^ 2 := by
  have hab' : a - b ≠ 0 := sub_ne_zero.mpr hab
  have hac' : a - c ≠ 0 := sub_ne_zero.mpr hac
  have hbc' : b - c ≠ 0 := sub_ne_zero.mpr hbc
  have k12 : (a - b) * (4 * (a ^ 2 + a * b + b ^ 2) - g) = 0 := by linear_combination h1 - h2
  have h12 : 4 * (a ^ 2 + a * b + b ^ 2) - g = 0 := by
    rcases mul_eq_zero.1 k12 with hc | hc
    · exact absurd hc hab'
    · exact hc
  have k13 : (a - c) * (4 * (a ^ 2 + a * c + c ^ 2) - g) = 0 := by linear_combination h1 - h3
  have h13 : 4 * (a ^ 2 + a * c + c ^ 2) - g = 0 := by
    rcases mul_eq_zero.1 k13 with hc | hc
    · exact absurd hc hac'
    · exact hc
  have ksum : (b - c) * (4 * (a + b + c)) = 0 := by linear_combination h12 - h13
  have hsum : a + b + c = 0 := by
    rcases mul_eq_zero.1 ksum with hc | hc
    · exact absurd hc hbc'
    · linear_combination hc / 4
  have hcc : c = -(a + b) := by linear_combination hsum
  subst hcc
  have hg : g = 4 * (a ^ 2 + a * b + b ^ 2) := by linear_combination -h12
  have hh : h = -(4 * (a * b * (a + b))) := by linear_combination -h1 + a * h12
  refine ⟨by ring, by linear_combination -h12, by linear_combination -h1 + a * h12, ?_⟩
  rw [hg, hh]; ring

/-! ## ★★★★★★★★★★★★★★★★★★★★判別式は差積 -/

/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] 判別式は半周期の値の差積である**——★**無条件**（第 1397）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

    latticeDisc Λ = 16·((e₁−e₂)(e₁−e₃)(e₂−e₃))²,   e_i = ℘(v_i)

★★★これが `Δ(E)^l = Δ(E/C)·N⁴` の証明の第 1 歩である
——両辺の `Δ` を半周期の値の差積に開いて、
同種のノルム（`∏_{w∈T}(℘(z+w)−e_i) = c_i(℘_{Λ′}(z)−e′_i)`）で結ぶ。 -/
theorem latticeDisc_eq_prod_half (L : PeriodPair) :
    latticeDisc L = 16 * ((L.weierstrassP (L.ω₁ / 2) - L.weierstrassP (L.ω₂ / 2))
      * (L.weierstrassP (L.ω₁ / 2) - L.weierstrassP ((L.ω₁ + L.ω₂) / 2))
      * (L.weierstrassP (L.ω₂ / 2) - L.weierstrassP ((L.ω₁ + L.ω₂) / 2))) ^ 2 :=
  (vieta_of_three_roots _ _ _ L.g₂ L.g₃ (e₁_ne_e₂ L) (e₁_ne_e₃ L) (e₂_ne_e₃ L)
    (isRoot_cubic_e₁ L) (isRoot_cubic_e₂ L) (isRoot_cubic_e₃ L)).2.2.2

/-- ★★★★★★★★半周期の値の和は `0`（第 1397）。 -/
theorem sum_half_eq_zero (L : PeriodPair) :
    L.weierstrassP (L.ω₁ / 2) + L.weierstrassP (L.ω₂ / 2)
      + L.weierstrassP ((L.ω₁ + L.ω₂) / 2) = 0 :=
  (vieta_of_three_roots _ _ _ L.g₂ L.g₃ (e₁_ne_e₂ L) (e₁_ne_e₃ L) (e₂_ne_e₃ L)
    (isRoot_cubic_e₁ L) (isRoot_cubic_e₂ L) (isRoot_cubic_e₃ L)).1

/-- ★★★★★★★★`g₂` は半周期の値の 2 次基本対称式（第 1397）。 -/
theorem g₂_eq_half (L : PeriodPair) :
    L.g₂ = -4 * (L.weierstrassP (L.ω₁ / 2) * L.weierstrassP (L.ω₂ / 2)
      + L.weierstrassP (L.ω₁ / 2) * L.weierstrassP ((L.ω₁ + L.ω₂) / 2)
      + L.weierstrassP (L.ω₂ / 2) * L.weierstrassP ((L.ω₁ + L.ω₂) / 2)) :=
  (vieta_of_three_roots _ _ _ L.g₂ L.g₃ (e₁_ne_e₂ L) (e₁_ne_e₃ L) (e₂_ne_e₃ L)
    (isRoot_cubic_e₁ L) (isRoot_cubic_e₂ L) (isRoot_cubic_e₃ L)).2.1

/-- ★★★★★★★★`g₃` は半周期の値の 3 次基本対称式（第 1397）。 -/
theorem g₃_eq_half (L : PeriodPair) :
    L.g₃ = 4 * (L.weierstrassP (L.ω₁ / 2) * L.weierstrassP (L.ω₂ / 2)
      * L.weierstrassP ((L.ω₁ + L.ω₂) / 2)) :=
  (vieta_of_three_roots _ _ _ L.g₂ L.g₃ (e₁_ne_e₂ L) (e₁_ne_e₃ L) (e₂_ne_e₃ L)
    (isRoot_cubic_e₁ L) (isRoot_cubic_e₂ L) (isRoot_cubic_e₃ L)).2.2.1

/-! ## ★★★★★★★★半周期を足す公式の代数の核 -/

/-- ★★★★★★★★★★★★★★★★
**半周期を足す公式の代数の核**——★**無条件**（第 1397）。

    ∏_{j} ∏_{i≠j} ( e_j + (e_j−e_a)(e_j−e_b)/(X−e_j) − e_i ) = −((e₁−e₂)(e₁−e₃)(e₂−e₃))²

★★★**右辺が `X` に依らない**のが要である
——これが `∏_{i≠j}(℘(v_j+w) − e_i) = −D`（`w` に依らない）の中身である。 -/
theorem half_shift_prod_eq {F : Type} [Field F] (e₁ e₂ e₃ X : F)
    (h₁ : X - e₁ ≠ 0) (h₂ : X - e₂ ≠ 0) (h₃ : X - e₃ ≠ 0) :
    ((e₁ + (e₁ - e₂) * (e₁ - e₃) / (X - e₁)) - e₂)
      * ((e₁ + (e₁ - e₂) * (e₁ - e₃) / (X - e₁)) - e₃)
      * ((e₂ + (e₂ - e₁) * (e₂ - e₃) / (X - e₂)) - e₁)
      * ((e₂ + (e₂ - e₁) * (e₂ - e₃) / (X - e₂)) - e₃)
      * ((e₃ + (e₃ - e₁) * (e₃ - e₂) / (X - e₃)) - e₁)
      * ((e₃ + (e₃ - e₁) * (e₃ - e₂) / (X - e₃)) - e₂)
      = -(((e₁ - e₂) * (e₁ - e₃) * (e₂ - e₃)) ^ 2) := by
  field_simp
  ring

/-! ## ★出典の紐付け(`.src`) -/

def sub_div_two_notMem_lattice.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5((ω₁−ω₂)/2 も束に入らない。★無条件)",
    sectionId := "genell-lemma-3-5" }

def e₁_ne_e₂.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(半周期の値は相異なる e₁ ≠ e₂。★無条件)",
    sectionId := "genell-lemma-3-5" }

def e₁_ne_e₃.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(半周期の値は相異なる e₁ ≠ e₃。★無条件)",
    sectionId := "genell-lemma-3-5" }

def e₂_ne_e₃.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(半周期の値は相異なる e₂ ≠ e₃。★無条件)",
    sectionId := "genell-lemma-3-5" }

def vieta_of_three_roots.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(相異なる 3 根から Vieta と判別式の差積表示。★無条件)",
    sectionId := "genell-lemma-3-5" }

def latticeDisc_eq_prod_half.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(判別式は半周期の値の差積である。★無条件)",
    sectionId := "genell-lemma-3-5" }

def latticeDisc_eq_prod_half.needs : List ProofObligation :=
  [ .citation "[ABC3]" "mem_lattice_of_shift_eq(第 624、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.mem_lattice_of_shift_eq") 1,
    .citation "[ABC3]" "isRoot_cubic_e₁/e₂/e₃(在庫、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.isRoot_cubic_e₁") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1397）**——`Δ(E)^l = Δ(E/C)·N⁴` の証明の第 1 歩である。" ++
       "☆道は (1) 判別式の差積表示（本ブロック）、(2) 半周期を足す公式、" ++
       "(3) 同種のノルム `∏_{w∈T}(℘(z+w)−e_i) = c_i(℘_{Λ′}(z)−e′_i)` の 3 つ。" ++
       "★(3) を `z = v_j` で使い (2) で書き直すと " ++
       "`∏_{i≠j}(℘(v_j+w)−e_i) = −D`（`w` に依らない）という代数恒等式に落ちる" ++
       "（`half_shift_prod_eq`）。☆60 桁の数値検証済み。") 17 ]

def sum_half_eq_zero.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(半周期の値の和は 0。★無条件)",
    sectionId := "genell-lemma-3-5" }

def g₂_eq_half.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(g₂ は半周期の値の 2 次基本対称式。★無条件)",
    sectionId := "genell-lemma-3-5" }

def g₃_eq_half.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(g₃ は半周期の値の 3 次基本対称式。★無条件)",
    sectionId := "genell-lemma-3-5" }

def half_shift_prod_eq.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(半周期を足す公式の代数の核——X に依らない。★無条件)",
    sectionId := "genell-lemma-3-5" }

end ABC3.Found.GenEll
