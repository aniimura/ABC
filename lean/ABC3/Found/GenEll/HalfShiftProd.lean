/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GenEll.HalfPeriodDisc
import ABC3.Meta.Claim

/-!
# 第 1398 ブロック —— **`∏_{i≠j}(℘(v_j+w) − e_i) = −D`（`w` に依らない）**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.17。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

## ★★★★★★★★★★★★★★★★★★★★これは何か——恒等式の証明の第 2 歩

第 1397 で `Δ = 16·((e₁−e₂)(e₁−e₃)(e₂−e₃))²`（`D` と書く）を得た。
★本ブロックは**半周期を足したときの `℘` の値**を計算し、

    ∏_{j} ∏_{i≠j} ( ℘(v_j + w) − e_i ) = −D     （`2w ∉ Λ` なら **`w` に依らない**）

を示す。☆`w` に依らないのが要である——最終形で `w ∈ T∖{0}` について積を取ると
`(−D)^{l−1} = D^{l−1}` になる。

## ★★★★★★★★道具は在庫だけ——**新しい解析は要らなかった**

| 段 | 道具 | 出どころ |
|---|---|---|
| `℘` の加法定理 | `weierstrassP_addition_general` | 在庫（第 655） |
| 半周期で `℘′ = 0` | `derivWeierstrassP_eq_zero_of_two_mem` | 在庫（第 605） |
| 単射性（`℘(v) ≠ ℘(w)` を出す） | `sub_or_add_mem_of_weierstrassP_eq` | 在庫（第 660） |
| `℘′² = 4∏(℘ − e_i)` | `derivWeierstrassP_sq` ＋ 第 1397 の Vieta | 在庫＋第 1397 |
| 最後の代数 | `field_simp; ring` | ★ |

★★★**半周期を足す公式**は加法定理から 3 行で出る:

    ℘(v + w) = ℘′(w)² / (4(℘(v) − ℘(w))²) − ℘(v) − ℘(w)     （`℘′(v) = 0` だから）

☆これに `℘′(w)² = 4(X−e₁)(X−e₂)(X−e₃)`（`X = ℘(w)`）を入れると

    ℘(v_j + w) = (X−e_a)(X−e_b)/(X−e_j) − e_j − X

となり、あとは**自由な有理式の恒等式**である（`e₁+e₂+e₃ = 0` だけ使う）。
-/

namespace ABC3.Found.GenEll

open PeriodPair ABC3.Meta

set_option maxRecDepth 8000

/-! ## ★★★★★★★★半周期を足したときの `℘` の値 -/

/-- ★★★★★★★★★★★★
**2-捩れ点を足したときの `℘` の値**——★**無条件**（第 1398）。

    ℘(v + w) = ℘′(w)² / (4(℘(v) − ℘(w))²) − ℘(v) − ℘(w)

☆加法定理に `℘′(v) = 0` を入れただけである。 -/
theorem weierstrassP_half_shift (P : PeriodPair) {v w : ℂ}
    (hv : v ∉ P.lattice) (h2v : 2 * v ∈ P.lattice) (h2w : 2 * w ∉ P.lattice) :
    P.weierstrassP (v + w)
      = (P.derivWeierstrassP w) ^ 2 / (4 * (P.weierstrassP v - P.weierstrassP w) ^ 2)
        - P.weierstrassP v - P.weierstrassP w := by
  have hw : w ∉ P.lattice := fun hc => h2w (by
    have h : (2 : ℂ) * w = w + w := by ring
    rw [h]; exact P.lattice.add_mem hc hc)
  have hvw : v + w ∉ P.lattice := by
    intro hc
    refine h2w ?_
    have he : (2 : ℂ) * w = ((v + w) + (v + w)) - 2 * v := by ring
    rw [he]
    exact P.lattice.sub_mem (P.lattice.add_mem hc hc) h2v
  have hne : P.weierstrassP v - P.weierstrassP w ≠ 0 := by
    intro hc
    have heq : P.weierstrassP v = P.weierstrassP w := by linear_combination hc
    rcases sub_or_add_mem_of_weierstrassP_eq P hv hw heq with hsub | hadd
    · refine h2w ?_
      have he : (2 : ℂ) * w = 2 * v - ((v - w) + (v - w)) := by ring
      rw [he]
      exact P.lattice.sub_mem h2v (P.lattice.add_mem hsub hsub)
    · exact hvw hadd
  have hdv : P.derivWeierstrassP v = 0 := derivWeierstrassP_eq_zero_of_two_mem P v h2v
  rw [weierstrassP_addition_general P hv hw hvw hne, hdv]
  field_simp
  ring

/-- ★★★★★★**2-捩れ点と非 2-捩れ点では `℘` の値が違う**（第 1398）。 -/
theorem weierstrassP_half_ne (P : PeriodPair) {v w : ℂ}
    (hv : v ∉ P.lattice) (h2v : 2 * v ∈ P.lattice) (h2w : 2 * w ∉ P.lattice) :
    P.weierstrassP w - P.weierstrassP v ≠ 0 := by
  have hw : w ∉ P.lattice := fun hc => h2w (by
    have h : (2 : ℂ) * w = w + w := by ring
    rw [h]; exact P.lattice.add_mem hc hc)
  have hvw : v + w ∉ P.lattice := by
    intro hc
    refine h2w ?_
    have he : (2 : ℂ) * w = ((v + w) + (v + w)) - 2 * v := by ring
    rw [he]
    exact P.lattice.sub_mem (P.lattice.add_mem hc hc) h2v
  intro hc
  have heq : P.weierstrassP v = P.weierstrassP w := by linear_combination -hc
  rcases sub_or_add_mem_of_weierstrassP_eq P hv hw heq with hsub | hadd
  · refine h2w ?_
    have he : (2 : ℂ) * w = 2 * v - ((v - w) + (v - w)) := by ring
    rw [he]
    exact P.lattice.sub_mem h2v (P.lattice.add_mem hsub hsub)
  · exact hvw hadd

/-- ★★★★★★★★**`℘′² = 4(X−e₁)(X−e₂)(X−e₃)`**——第 1397 の Vieta から（第 1398）。 -/
theorem derivSq_factor (P : PeriodPair) {w : ℂ} (hw : w ∉ P.lattice) :
    P.derivWeierstrassP w ^ 2
      = 4 * (P.weierstrassP w - P.weierstrassP (P.ω₁ / 2))
        * (P.weierstrassP w - P.weierstrassP (P.ω₂ / 2))
        * (P.weierstrassP w - P.weierstrassP ((P.ω₁ + P.ω₂) / 2)) := by
  rw [P.derivWeierstrassP_sq w hw, g₂_eq_half P, g₃_eq_half P]
  linear_combination (4 * P.weierstrassP w ^ 2) * (sum_half_eq_zero P)

/-! ## ★★★★代数の橋 -/

theorem shift_val_eq₁ (a b c X : ℂ) (h : X - a ≠ 0) :
    4 * (X - a) * (X - b) * (X - c) / (4 * (a - X) ^ 2) - a - X
      = (X - b) * (X - c) / (X - a) - a - X := by
  rw [show (a - X) ^ 2 = (X - a) ^ 2 by ring]; field_simp

theorem shift_val_eq₂ (a b c X : ℂ) (h : X - b ≠ 0) :
    4 * (X - a) * (X - b) * (X - c) / (4 * (b - X) ^ 2) - b - X
      = (X - a) * (X - c) / (X - b) - b - X := by
  rw [show (b - X) ^ 2 = (X - b) ^ 2 by ring]; field_simp

theorem shift_val_eq₃ (a b c X : ℂ) (h : X - c ≠ 0) :
    4 * (X - a) * (X - b) * (X - c) / (4 * (c - X) ^ 2) - c - X
      = (X - a) * (X - b) / (X - c) - c - X := by
  rw [show (c - X) ^ 2 = (X - c) ^ 2 by ring]; field_simp

/-- ★★★★★★★★★★★★★★★★
**代数の核**——★**無条件**（第 1398）。

☆`e₁+e₂+e₃ = 0` だけを使う自由な有理式の恒等式であり、
**右辺が `X` に依らない**のが要である。 -/
theorem prod_half_shift_alg (a b c X : ℂ)
    (hsum : a + b + c = 0) (ha : X - a ≠ 0) (hb : X - b ≠ 0) (hc : X - c ≠ 0) :
    (((X - b) * (X - c) / (X - a) - a - X) - b)
      * (((X - b) * (X - c) / (X - a) - a - X) - c)
      * (((X - a) * (X - c) / (X - b) - b - X) - a)
      * (((X - a) * (X - c) / (X - b) - b - X) - c)
      * (((X - a) * (X - b) / (X - c) - c - X) - a)
      * (((X - a) * (X - b) / (X - c) - c - X) - b)
      = -(((a - b) * (a - c) * (b - c)) ^ 2) := by
  have hc' : c = -(a + b) := by linear_combination hsum
  subst hc'
  field_simp
  ring

/-! ## ★★★★★★★★★★★★★★★★★★★★結論 -/

/-- ★★★★★★★★★★★★★★★★★★★★
**[GenEll] `∏_{i≠j}(℘(v_j+w) − e_i) = −D`**——★**`2w ∉ Λ` だけ**（第 1398）。

原文 (GenEll p.17):
> Lemma 3.5. (Global Rank One Subgroups of l-Torsion) Let

★★★**右辺は `w` に依らない**——これが `Δ(E)^l = Δ(E/C)·N⁴` の証明の心臓である。
☆最終形では `w ∈ T∖{0}`（`l−1` 個）について積を取り `(−D)^{l−1} = D^{l−1}` を使う。 -/
theorem prod_half_shift (P : PeriodPair) {w : ℂ} (h2w : 2 * w ∉ P.lattice) :
    (P.weierstrassP (P.ω₁ / 2 + w) - P.weierstrassP (P.ω₂ / 2))
      * (P.weierstrassP (P.ω₁ / 2 + w) - P.weierstrassP ((P.ω₁ + P.ω₂) / 2))
      * (P.weierstrassP (P.ω₂ / 2 + w) - P.weierstrassP (P.ω₁ / 2))
      * (P.weierstrassP (P.ω₂ / 2 + w) - P.weierstrassP ((P.ω₁ + P.ω₂) / 2))
      * (P.weierstrassP ((P.ω₁ + P.ω₂) / 2 + w) - P.weierstrassP (P.ω₁ / 2))
      * (P.weierstrassP ((P.ω₁ + P.ω₂) / 2 + w) - P.weierstrassP (P.ω₂ / 2))
      = -(((P.weierstrassP (P.ω₁ / 2) - P.weierstrassP (P.ω₂ / 2))
          * (P.weierstrassP (P.ω₁ / 2) - P.weierstrassP ((P.ω₁ + P.ω₂) / 2))
          * (P.weierstrassP (P.ω₂ / 2) - P.weierstrassP ((P.ω₁ + P.ω₂) / 2))) ^ 2) := by
  have hw : w ∉ P.lattice := fun hc => h2w (by
    have h : (2 : ℂ) * w = w + w := by ring
    rw [h]; exact P.lattice.add_mem hc hc)
  have h2v₁ : 2 * (P.ω₁ / 2) ∈ P.lattice := by
    rw [show 2 * (P.ω₁ / 2) = P.ω₁ by ring]; exact P.ω₁_mem_lattice
  have h2v₂ : 2 * (P.ω₂ / 2) ∈ P.lattice := by
    rw [show 2 * (P.ω₂ / 2) = P.ω₂ by ring]; exact P.ω₂_mem_lattice
  have h2v₃ : 2 * ((P.ω₁ + P.ω₂) / 2) ∈ P.lattice := by
    rw [show 2 * ((P.ω₁ + P.ω₂) / 2) = P.ω₁ + P.ω₂ by ring]
    exact P.lattice.add_mem P.ω₁_mem_lattice P.ω₂_mem_lattice
  have hne₁ := weierstrassP_half_ne P P.ω₁_div_two_notMem_lattice h2v₁ h2w
  have hne₂ := weierstrassP_half_ne P P.ω₂_div_two_notMem_lattice h2v₂ h2w
  have hne₃ := weierstrassP_half_ne P (add_div_two_notMem_lattice P) h2v₃ h2w
  have h1 : P.weierstrassP (P.ω₁ / 2 + w)
      = (P.weierstrassP w - P.weierstrassP (P.ω₂ / 2))
        * (P.weierstrassP w - P.weierstrassP ((P.ω₁ + P.ω₂) / 2))
        / (P.weierstrassP w - P.weierstrassP (P.ω₁ / 2))
        - P.weierstrassP (P.ω₁ / 2) - P.weierstrassP w := by
    rw [weierstrassP_half_shift P P.ω₁_div_two_notMem_lattice h2v₁ h2w, derivSq_factor P hw]
    exact shift_val_eq₁ _ _ _ _ hne₁
  have h2 : P.weierstrassP (P.ω₂ / 2 + w)
      = (P.weierstrassP w - P.weierstrassP (P.ω₁ / 2))
        * (P.weierstrassP w - P.weierstrassP ((P.ω₁ + P.ω₂) / 2))
        / (P.weierstrassP w - P.weierstrassP (P.ω₂ / 2))
        - P.weierstrassP (P.ω₂ / 2) - P.weierstrassP w := by
    rw [weierstrassP_half_shift P P.ω₂_div_two_notMem_lattice h2v₂ h2w, derivSq_factor P hw]
    exact shift_val_eq₂ _ _ _ _ hne₂
  have h3 : P.weierstrassP ((P.ω₁ + P.ω₂) / 2 + w)
      = (P.weierstrassP w - P.weierstrassP (P.ω₁ / 2))
        * (P.weierstrassP w - P.weierstrassP (P.ω₂ / 2))
        / (P.weierstrassP w - P.weierstrassP ((P.ω₁ + P.ω₂) / 2))
        - P.weierstrassP ((P.ω₁ + P.ω₂) / 2) - P.weierstrassP w := by
    rw [weierstrassP_half_shift P (add_div_two_notMem_lattice P) h2v₃ h2w,
      derivSq_factor P hw]
    exact shift_val_eq₃ _ _ _ _ hne₃
  rw [h1, h2, h3]
  exact prod_half_shift_alg _ _ _ _ (sum_half_eq_zero P) hne₁ hne₂ hne₃

/-! ## ★出典の紐付け(`.src`) -/

def weierstrassP_half_shift.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2-捩れ点を足したときの ℘ の値。★無条件)",
    sectionId := "genell-lemma-3-5" }

def weierstrassP_half_ne.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(2-捩れ点と非 2-捩れ点では ℘ の値が違う。★無条件)",
    sectionId := "genell-lemma-3-5" }

def derivSq_factor.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(℘′² = 4(X−e₁)(X−e₂)(X−e₃)。★無条件)",
    sectionId := "genell-lemma-3-5" }

def shift_val_eq₁.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(半周期を足す公式の代数の橋 その 1。★無条件)",
    sectionId := "genell-lemma-3-5" }

def shift_val_eq₂.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(半周期を足す公式の代数の橋 その 2。★無条件)",
    sectionId := "genell-lemma-3-5" }

def shift_val_eq₃.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(半周期を足す公式の代数の橋 その 3。★無条件)",
    sectionId := "genell-lemma-3-5" }

def prod_half_shift_alg.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(代数の核——右辺が X に依らない。★無条件)",
    sectionId := "genell-lemma-3-5" }

def prod_half_shift.src : Source :=
  { paper := "GenEll", pdfPage := 17,
    item := "Lemma 3.5(∏_{i≠j}(℘(v_j+w) − e_i) = −D——w に依らない。★2w ∉ Λ)",
    sectionId := "genell-lemma-3-5" }

def prod_half_shift.needs : List ProofObligation :=
  [ .citation "[ABC3]" "weierstrassP_addition_general(第 655、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.weierstrassP_addition_general") 1,
    .citation "[ABC3]" "latticeDisc_eq_prod_half・sum_half_eq_zero(第 1397、証明済み)"
      (.inProject "ABC3" "ABC3.Found.GenEll.latticeDisc_eq_prod_half") 1,
    .implicitStep
      ("★★★★**2026-09-02（第 1398）**——`Δ(E)^l = Δ(E/C)·N⁴` の証明の第 2 歩である。" ++
       "☆**新しい解析は要らなかった**——加法定理（第 655）と単射性（第 660）と" ++
       "第 1397 の Vieta だけで出た。" ++
       "★★★残るのは**同種のノルム** " ++
       "`∏_{w∈T}(℘(z+w) − e_i) = c_i·(℘_{Λ′}(z) − e′_i)` の 1 本だけである" ++
       "——`elliptic_liouville`（第 598、在庫）と極の解析で出る見込み。") 17 ]

end ABC3.Found.GenEll
