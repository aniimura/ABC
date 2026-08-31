/-
Copyright (c) 2026 ABC3 Project. All rights reserved.
-/
import ABC3.Found.GaloisRep.MuCharSum
import Mathlib.Tactic.LinearCombination

/-!
# Galois (G6) 第 847 ブロック —— **★★★★★★★★★★`∑_ζ ζ^a/(1−ζ)^k` の機械**（`Found`）

原典: S. Mochizuki, *Arithmetic Elliptic Curves in General Position* [GenEll]、物理 p.15。

原文 (GenEll p.15):
> parameter qE of E satisfies the relation qE = qEl ; in particular, we have

## ★★★★★★★★★★★★★★★★★★★★これは何か

`sum_mu_dyterm`（第 846 の葉）は `240·∑_{ζ≠1} ζ²(2+ζ)/(1−ζ)⁴ = l⁴ − 1` であり、
`p₄ = ∑_{ζ≠1} 1/(1−ζ)⁴` が要る。★これまでの道具（`sum_mu_frac` 等）は
`p₁, p₂, p₃` までで、`p₄` には新しい種が要る。

★★**種は要らない**——次の 3 つの性質だけで**すべての `k` が決まる**:

    `A^{(k)}_a := ∑_{ζ≠1} ζ^a/(1−ζ)^k`

| 性質 | 内容 | 根拠 |
|---|---|---|
| (i) | `A^{(k+1)}_a − A^{(k+1)}_{a+1} = A^{(k)}_a` | `ζ^a·(1−ζ)/(1−ζ)^{k+1} = ζ^a/(1−ζ)^k` |
| (ii) | `∑_{a<l} A^{(k)}_a = 0` | ★★**`∑_{a<l} ζ^a = 0`**（`ζ ≠ 1`） |
| (iii) | `A^{(0)}_a = l·[l∣a] − 1` | 指標和 |

(i) を繰り返すと `A^{(k+1)}_a = A^{(k+1)}_0 − ∑_{j<a} A^{(k)}_j`。これを (ii) に入れると

    ★★★★★★`l · A^{(k+1)}_0 = ∑_{j<l} (l − 1 − j) · A^{(k)}_j`

——★**下の層だけから上の層の定数項が決まる漸化式**である。

## ★★★★★★★★★★★★★★★★★★★★★★★★重み付きの和で登る（次のブロックの計画）

一般の `a` での閉じた式を作ると `∑ j^m` が増えていく。★そこで**重み付きの和**

    `M_m^{(k)} := ∑_{i<l} C(l−1−i, m) · A^{(k)}_i`

を使う。ホッケースティック `∑_{t<l} C(t,m) = C(l,m+1)` と (i)(ii) から

| 式 | 内容 |
|---|---|
| `M_0^{(k)} = 0` | (ii) そのもの |
| ★★`M_m^{(k+1)} = C(l,m+1)·A^{(k+1)}_0 − M_{m+1}^{(k)}` | (i) を代入してホッケースティック |
| `M_m^{(0)} = l·C(l−1,m) − C(l,m+1)` | (iii) から |
| ★`l·A^{(k+1)}_0 = M_1^{(k)}` | `M_0^{(k+1)} = 0` の別の書き方 |

★★★**一般の `a` の式は一切要らない**。`p₄` までの鎖は有限で短い:

```
M_4^{(0)} → M_3^{(1)} → M_2^{(2)} → M_1^{(3)} → A^{(4)}_0
M_3^{(0)} → M_2^{(1)} → M_1^{(2)} → A^{(3)}_0
M_2^{(0)} → M_1^{(1)} → A^{(2)}_0
M_1^{(0)} → A^{(1)}_0
```

☆検算: `A^{(1)}_0 = [l·C(l−1,1) − C(l,2)]/l = (l−1) − (l−1)/2 = (l−1)/2` ✓
☆`A^{(2)}_0 = [C(l,2)(l−1)/2 − l·C(l−1,2) + C(l,3)]/l = (l−1)(5−l)/12` ✓

## ★★★検算（`l = 2, 3, 5, 7, 11, 13, 17` で数値確認、2026-08-31）

    `p₁ = (l−1)/2`
    `p₂ = (l−1)(5−l)/12`
    `p₃ = (l−1)(3−l)/8`
    `p₄ = (l−1)(l³ + l² − 109l + 251)/720`

そして `3p₄ − 7p₃ + 5p₂ − p₁ = (l⁴−1)/240` が厳密に成り立つ
（`x = 1/(1−ζ)` と置くと `ζ²(2+ζ)/(1−ζ)⁴ = 3x⁴ − 7x³ + 5x² − x`）。

## ★本ブロックで取れるもの

| 定義・定理 | 内容 |
|---|---|
| `muPow` | `A^{(k)}_a = ∑_{ζ≠1} ζ^a/(1−ζ)^k` |
| `muPow_succ_sub` | ★(i) |
| `muPow_sum_range` | ★★(ii)——`∑_{a<l} ζ^a = 0` から |
| `muPow_zero_level` | (iii) |
| `muPow_shift` | `A^{(k+1)}_a = A^{(k+1)}_0 − ∑_{j<a}A^{(k)}_j` |
| `muPow_zero_rec` | ★★★★★★**漸化式** |
| `sum_range_sum_range_lt` | `∑_{a<n}∑_{j<a}f j = ∑_{j<n}(n−1−j)f j` |
-/

namespace ABC3.Found.GaloisRep

open Finset

variable {F : Type*} [Field F]

/-! ## ★★★★★機械の定義 -/

/-- ★★★★★★**`A^{(k)}_a = ∑_{ζ≠1} ζ^a/(1−ζ)^k`**。 -/
noncomputable def muPow (l : ℕ) (ζ : F) (k a : ℕ) : F :=
  ∑ i ∈ (range l).erase 0, (ζ ^ i) ^ a * ((1 - ζ ^ i)⁻¹) ^ k

theorem one_sub_zeta_ne_zero {l : ℕ} {ζ : F} (hζ : IsPrimitiveRoot ζ l)
    {i : ℕ} (hi : i ∈ (range l).erase 0) : (1 : F) - ζ ^ i ≠ 0 := by
  have hi0 : i ≠ 0 := (Finset.mem_erase.1 hi).1
  have hil : i < l := Finset.mem_range.1 (Finset.mem_erase.1 hi).2
  intro h
  exact hζ.pow_ne_one_of_pos_of_lt hi0 hil (sub_eq_zero.1 h).symm

/-- ★★★**(i) 一段下げる**——`ζ^a·(1−ζ)/(1−ζ)^{k+1} = ζ^a/(1−ζ)^k`。 -/
theorem muPow_succ_sub {l : ℕ} {ζ : F} (hζ : IsPrimitiveRoot ζ l) (k a : ℕ) :
    muPow l ζ (k + 1) a - muPow l ζ (k + 1) (a + 1) = muPow l ζ k a := by
  simp only [muPow, ← Finset.sum_sub_distrib]
  refine Finset.sum_congr rfl (fun i hi => ?_)
  have hinv : (1 - ζ ^ i) * ((1 - ζ ^ i)⁻¹) = 1 :=
    mul_inv_cancel₀ (one_sub_zeta_ne_zero hζ hi)
  linear_combination ((ζ ^ i) ^ a * ((1 - ζ ^ i)⁻¹) ^ k) * hinv

/-- ★★★★★★**(ii) `a` について一周すると消える**——`∑_{a<l} ζ^a = 0` から。 -/
theorem muPow_sum_range {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) (k : ℕ) :
    ∑ a ∈ range l, muPow l ζ k a = 0 := by
  simp only [muPow]
  rw [Finset.sum_comm]
  refine Finset.sum_eq_zero (fun i hi => ?_)
  obtain ⟨-, hsum⟩ := zeta_pow_facts hl hζ hi
  rw [← Finset.sum_mul, hsum, zero_mul]

/-- ★★(iii) 一番下の層。 -/
theorem muPow_zero_level {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) (a : ℕ) :
    muPow l ζ 0 a = (if l ∣ a then (l : F) else 0) - 1 := by
  simp only [muPow, pow_zero, mul_one]
  rw [← sum_mu_pow_erase_zero hl hζ a]
  exact Finset.sum_congr rfl (fun i _ => by rw [← pow_mul])

/-- ★★★★**(i) を繰り返した形**。 -/
theorem muPow_shift {l : ℕ} {ζ : F} (hζ : IsPrimitiveRoot ζ l) (k a : ℕ) :
    muPow l ζ (k + 1) a = muPow l ζ (k + 1) 0 - ∑ j ∈ range a, muPow l ζ k j := by
  induction a with
  | zero => simp
  | succ n ih =>
      have h := muPow_succ_sub hζ k n
      rw [Finset.sum_range_succ]
      linear_combination -h + ih

/-- ★`∑_{a<n}∑_{j<a} f j = ∑_{j<n} (n − 1 − j)·f j`。 -/
theorem sum_range_sum_range_lt (n : ℕ) (f : ℕ → F) :
    ∑ a ∈ range n, ∑ j ∈ range a, f j = ∑ j ∈ range n, ((n : F) - 1 - (j : F)) * f j := by
  induction n with
  | zero => simp
  | succ m ih =>
      rw [Finset.sum_range_succ, ih, Finset.sum_range_succ]
      push_cast
      rw [show ((m : F) + 1 - 1 - (m : F)) = 0 by ring, zero_mul, add_zero,
        ← Finset.sum_add_distrib]
      refine Finset.sum_congr rfl (fun j _ => ?_)
      ring

/-- ★★★★★★★★★★**漸化式**——下の層だけから上の層の定数項が決まる。

    `l · A^{(k+1)}_0 = ∑_{j<l} (l − 1 − j) · A^{(k)}_j` -/
theorem muPow_zero_rec {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) (k : ℕ) :
    (l : F) * muPow l ζ (k + 1) 0
      = ∑ j ∈ range l, ((l : F) - 1 - (j : F)) * muPow l ζ k j := by
  have h0 := muPow_sum_range hl hζ (k + 1)
  rw [Finset.sum_congr rfl (fun a _ => muPow_shift hζ k a)] at h0
  rw [Finset.sum_sub_distrib, Finset.sum_const, Finset.card_range, nsmul_eq_mul,
    sum_range_sum_range_lt] at h0
  linear_combination h0

/-- ★`i ↦ l − i` は `(range l).erase 0` の全単射である。 -/
theorem sum_erase_reflect {M : Type*} [AddCommMonoid M] {l : ℕ} (f : ℕ → M) :
    ∑ i ∈ (range l).erase 0, f (l - i) = ∑ i ∈ (range l).erase 0, f i := by
  refine Finset.sum_nbij' (fun i => l - i) (fun i => l - i) ?_ ?_ ?_ ?_ ?_
  · intro a ha
    have ha0 : a ≠ 0 := (Finset.mem_erase.1 ha).1
    have hal : a < l := Finset.mem_range.1 (Finset.mem_erase.1 ha).2
    exact Finset.mem_erase.2 ⟨by omega, Finset.mem_range.2 (by omega)⟩
  · intro a ha
    have ha0 : a ≠ 0 := (Finset.mem_erase.1 ha).1
    have hal : a < l := Finset.mem_range.1 (Finset.mem_erase.1 ha).2
    exact Finset.mem_erase.2 ⟨by omega, Finset.mem_range.2 (by omega)⟩
  · intro a ha
    have ha0 : a ≠ 0 := (Finset.mem_erase.1 ha).1
    have hal : a < l := Finset.mem_range.1 (Finset.mem_erase.1 ha).2
    omega
  · intro a ha
    have ha0 : a ≠ 0 := (Finset.mem_erase.1 ha).1
    have hal : a < l := Finset.mem_range.1 (Finset.mem_erase.1 ha).2
    omega
  · intro a _
    rfl

/-! ## ★★★★★★★★重み付きの和 `M_m^{(k)}` -/

/-- ★ホッケースティック `∑_{t<l} C(t,m) = C(l,m+1)`。 -/
theorem sum_range_choose_eq (l m : ℕ) : ∑ t ∈ range l, t.choose m = l.choose (m + 1) := by
  induction l with
  | zero => simp
  | succ n ih =>
      rw [Finset.sum_range_succ, ih, Nat.choose_succ_succ]
      simp only [Nat.succ_eq_add_one]
      omega

/-- ★反転したホッケースティック。 -/
theorem sum_range_choose_reflect (l m : ℕ) :
    ∑ i ∈ range l, ((Nat.choose (l - 1 - i) m : ℕ) : F) = ((l.choose (m + 1) : ℕ) : F) := by
  rw [Finset.sum_range_reflect (fun t => ((Nat.choose t m : ℕ) : F)) l]
  rw [← Nat.cast_sum, sum_range_choose_eq]

/-- ★重み付きの二重和の入れ替え。 -/
theorem sum_range_weighted_sum (n : ℕ) (w f : ℕ → F) :
    ∑ i ∈ range n, w i * ∑ j ∈ range i, f j
      = ∑ j ∈ range n, f j * ∑ i ∈ Finset.Ico (j + 1) n, w i := by
  induction n with
  | zero => simp
  | succ m ih =>
      rw [Finset.sum_range_succ, ih, Finset.sum_range_succ]
      have hlast : ∑ i ∈ Finset.Ico (m + 1) (m + 1), w i = 0 := by simp
      rw [hlast, mul_zero, add_zero]
      have hstep : ∀ j ∈ range m, f j * ∑ i ∈ Finset.Ico (j + 1) (m + 1), w i
          = f j * ∑ i ∈ Finset.Ico (j + 1) m, w i + w m * f j := by
        intro j hj
        have hjm : j + 1 ≤ m := Finset.mem_range.1 hj
        rw [Finset.sum_Ico_succ_top hjm]
        ring
      rw [Finset.sum_congr rfl hstep, Finset.sum_add_distrib, ← Finset.mul_sum]

/-- ★`Ico` 上のホッケースティック。 -/
theorem ico_choose_sum (l m j : ℕ) :
    ∑ i ∈ Finset.Ico (j + 1) l, ((Nat.choose (l - 1 - i) m : ℕ) : F)
      = ((Nat.choose (l - 1 - j) (m + 1) : ℕ) : F) := by
  have hn : l - (j + 1) = l - 1 - j := by omega
  rw [Finset.sum_Ico_eq_sum_range, hn]
  have hterm : ∀ i ∈ range (l - 1 - j),
      ((Nat.choose (l - 1 - (j + 1 + i)) m : ℕ) : F)
        = ((Nat.choose ((l - 1 - j) - 1 - i) m : ℕ) : F) := by
    intro i _
    congr 2
    omega
  rw [Finset.sum_congr rfl hterm]
  exact sum_range_choose_reflect (l - 1 - j) m

/-- ★★★★★★**重み付きの和** `M_m^{(k)} = ∑_{i<l} C(l−1−i, m)·A^{(k)}_i`。 -/
noncomputable def muM (l : ℕ) (ζ : F) (k m : ℕ) : F :=
  ∑ i ∈ range l, ((Nat.choose (l - 1 - i) m : ℕ) : F) * muPow l ζ k i

/-- ★★**`M_0^{(k)} = 0`**——(ii) そのもの。 -/
theorem muM_weight_zero {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) (k : ℕ) :
    muM l ζ k 0 = 0 := by
  simp only [muM, Nat.choose_zero_right, Nat.cast_one, one_mul]
  exact muPow_sum_range hl hζ k

/-- ★★★★★★★★**層をひとつ上げる**:
`M_m^{(k+1)} = C(l,m+1)·A^{(k+1)}_0 − M_{m+1}^{(k)}`。 -/
theorem muM_rec {l : ℕ} {ζ : F} (hζ : IsPrimitiveRoot ζ l) (k m : ℕ) :
    muM l ζ (k + 1) m
      = ((l.choose (m + 1) : ℕ) : F) * muPow l ζ (k + 1) 0 - muM l ζ k (m + 1) := by
  have hexp : ∀ i ∈ range l,
      ((Nat.choose (l - 1 - i) m : ℕ) : F) * muPow l ζ (k + 1) i
        = ((Nat.choose (l - 1 - i) m : ℕ) : F) * muPow l ζ (k + 1) 0
          - ((Nat.choose (l - 1 - i) m : ℕ) : F) * ∑ j ∈ range i, muPow l ζ k j := by
    intro i _
    rw [muPow_shift hζ k i]
    ring
  simp only [muM]
  rw [Finset.sum_congr rfl hexp, Finset.sum_sub_distrib, ← Finset.sum_mul,
    sum_range_choose_reflect, sum_range_weighted_sum]
  congr 1
  refine Finset.sum_congr rfl (fun j _ => ?_)
  rw [ico_choose_sum l m j]
  ring

/-- ★★★**一番下の層**: `M_m^{(0)} = l·C(l−1,m) − C(l,m+1)`。 -/
theorem muM_level_zero {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) (m : ℕ) :
    muM l ζ 0 m
      = (l : F) * ((Nat.choose (l - 1) m : ℕ) : F) - ((l.choose (m + 1) : ℕ) : F) := by
  have hif : ∀ i ∈ range l,
      ((Nat.choose (l - 1 - i) m : ℕ) : F) * muPow l ζ 0 i
        = (if i = 0 then (l : F) * ((Nat.choose (l - 1) m : ℕ) : F) else 0)
          - ((Nat.choose (l - 1 - i) m : ℕ) : F) := by
    intro i hi
    have hil : i < l := Finset.mem_range.1 hi
    rw [muPow_zero_level hl hζ i]
    by_cases h0 : i = 0
    · subst h0
      simp
      ring
    · rw [if_neg (fun h : l ∣ i => h0 (Nat.eq_zero_of_dvd_of_lt h hil)), if_neg h0]
      ring
  simp only [muM]
  rw [Finset.sum_congr rfl hif, Finset.sum_sub_distrib, sum_range_choose_reflect]
  congr 1
  rw [Finset.sum_ite_eq' (range l) 0 (fun _ => (l : F) * ((Nat.choose (l - 1) m : ℕ) : F))]
  rw [if_pos (Finset.mem_range.2 hl.pos)]

/-- ★★★★★★**定数項は `M_1` から出る**: `l·A^{(k+1)}_0 = M_1^{(k)}`。 -/
theorem muPow_zero_from_M {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) (k : ℕ) :
    (l : F) * muPow l ζ (k + 1) 0 = muM l ζ k 1 := by
  have h := muM_rec (l := l) hζ k 0
  rw [muM_weight_zero hl hζ (k + 1), Nat.choose_one_right] at h
  linear_combination -h

/-! ## ★★★★二項係数を多項式として見る -/

theorem cast_choose_two' (n : ℕ) :
    2 * ((n.choose 2 : ℕ) : F) = (n : F) * ((n : F) - 1) := by
  induction n with
  | zero => simp
  | succ k ih =>
      rw [Nat.choose_succ_succ, Nat.choose_one_right]
      push_cast
      push_cast at ih
      linear_combination ih

theorem cast_choose_three' (n : ℕ) :
    6 * ((n.choose 3 : ℕ) : F) = (n : F) * ((n : F) - 1) * ((n : F) - 2) := by
  induction n with
  | zero => simp
  | succ k ih =>
      rw [Nat.choose_succ_succ]
      push_cast
      push_cast at ih
      linear_combination ih + 3 * cast_choose_two' (F := F) k

theorem cast_choose_four' (n : ℕ) :
    24 * ((n.choose 4 : ℕ) : F)
      = (n : F) * ((n : F) - 1) * ((n : F) - 2) * ((n : F) - 3) := by
  induction n with
  | zero => simp
  | succ k ih =>
      rw [Nat.choose_succ_succ]
      push_cast
      push_cast at ih
      linear_combination ih + 4 * cast_choose_three' (F := F) k

theorem cast_choose_five' (n : ℕ) :
    120 * ((n.choose 5 : ℕ) : F)
      = (n : F) * ((n : F) - 1) * ((n : F) - 2) * ((n : F) - 3) * ((n : F) - 4) := by
  induction n with
  | zero => simp
  | succ k ih =>
      rw [Nat.choose_succ_succ]
      push_cast
      push_cast at ih
      linear_combination ih + 5 * cast_choose_four' (F := F) k

section CharZeroChain

variable [CharZero F]

theorem cast_choose_two_eq (n : ℕ) :
    ((n.choose 2 : ℕ) : F) = (n : F) * ((n : F) - 1) / 2 := by
  linear_combination cast_choose_two' (F := F) n / 2

theorem cast_choose_three_eq (n : ℕ) :
    ((n.choose 3 : ℕ) : F) = (n : F) * ((n : F) - 1) * ((n : F) - 2) / 6 := by
  linear_combination cast_choose_three' (F := F) n / 6

theorem cast_choose_four_eq (n : ℕ) :
    ((n.choose 4 : ℕ) : F)
      = (n : F) * ((n : F) - 1) * ((n : F) - 2) * ((n : F) - 3) / 24 := by
  linear_combination cast_choose_four' (F := F) n / 24

theorem cast_choose_five_eq (n : ℕ) :
    ((n.choose 5 : ℕ) : F)
      = (n : F) * ((n : F) - 1) * ((n : F) - 2) * ((n : F) - 3) * ((n : F) - 4) / 120 := by
  linear_combination cast_choose_five' (F := F) n / 120

/-! ## ★★★★★★★★★★`p₁, p₂, p₃, p₄` を登る -/

theorem cast_pred {l : ℕ} (hl : l.Prime) : ((l - 1 : ℕ) : F) = (l : F) - 1 := by
  have h1 : (1 : ℕ) ≤ l := hl.one_lt.le
  rw [Nat.cast_sub h1, Nat.cast_one]

theorem cast_prime_ne_zero {l : ℕ} (hl : l.Prime) : ((l : F)) ≠ 0 :=
  Nat.cast_ne_zero.2 hl.ne_zero

/-- `M_1^{(0)} = l(l−1)/2`。 -/
theorem muM_zero_one {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    muM l ζ 0 1 = (l : F) * ((l : F) - 1) / 2 := by
  rw [muM_level_zero hl hζ 1, Nat.choose_one_right, cast_pred (F := F) hl,
    cast_choose_two_eq (F := F) l]
  ring

/-- `M_2^{(0)} = l(l−1)(l−2)/3`。 -/
theorem muM_zero_two {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    muM l ζ 0 2 = (l : F) * ((l : F) - 1) * ((l : F) - 2) / 3 := by
  rw [muM_level_zero hl hζ 2, cast_choose_two_eq (F := F) (l - 1),
    cast_choose_three_eq (F := F) l, cast_pred (F := F) hl]
  ring

/-- `M_3^{(0)} = l(l−1)(l−2)(l−3)/8`。 -/
theorem muM_zero_three {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    muM l ζ 0 3 = (l : F) * ((l : F) - 1) * ((l : F) - 2) * ((l : F) - 3) / 8 := by
  rw [muM_level_zero hl hζ 3, cast_choose_three_eq (F := F) (l - 1),
    cast_choose_four_eq (F := F) l, cast_pred (F := F) hl]
  ring

/-- `M_4^{(0)} = l(l−1)(l−2)(l−3)(l−4)/30`。 -/
theorem muM_zero_four {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    muM l ζ 0 4
      = (l : F) * ((l : F) - 1) * ((l : F) - 2) * ((l : F) - 3) * ((l : F) - 4) / 30 := by
  rw [muM_level_zero hl hζ 4, cast_choose_four_eq (F := F) (l - 1),
    cast_choose_five_eq (F := F) l, cast_pred (F := F) hl]
  ring

/-- ★★★★**`p₁ = (l−1)/2`**。 -/
theorem muPow_one_zero {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    muPow l ζ 1 0 = ((l : F) - 1) / 2 := by
  refine mul_left_cancel₀ (cast_prime_ne_zero (F := F) hl) ?_
  rw [muPow_zero_from_M hl hζ 0, muM_zero_one hl hζ]
  ring

/-- `M_1^{(1)} = l(l−1)(5−l)/12`。 -/
theorem muM_one_one {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    muM l ζ 1 1 = (l : F) * ((l : F) - 1) * (5 - (l : F)) / 12 := by
  rw [muM_rec hζ 0 1, muPow_one_zero hl hζ, muM_zero_two hl hζ,
    cast_choose_two_eq (F := F) l]
  ring

/-- ★★★★**`p₂ = (l−1)(5−l)/12`**。 -/
theorem muPow_two_zero {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    muPow l ζ 2 0 = ((l : F) - 1) * (5 - (l : F)) / 12 := by
  refine mul_left_cancel₀ (cast_prime_ne_zero (F := F) hl) ?_
  rw [muPow_zero_from_M hl hζ 1, muM_one_one hl hζ]
  ring

/-- `M_2^{(1)} = l(l−1)(l−2)(7−l)/24`。 -/
theorem muM_one_two {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    muM l ζ 1 2 = (l : F) * ((l : F) - 1) * ((l : F) - 2) * (7 - (l : F)) / 24 := by
  rw [muM_rec hζ 0 2, muPow_one_zero hl hζ, muM_zero_three hl hζ,
    cast_choose_three_eq (F := F) l]
  ring

/-- `M_1^{(2)} = l(l−1)(3−l)/8`。 -/
theorem muM_two_one {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    muM l ζ 2 1 = (l : F) * ((l : F) - 1) * (3 - (l : F)) / 8 := by
  rw [muM_rec hζ 1 1, muPow_two_zero hl hζ, muM_one_two hl hζ,
    cast_choose_two_eq (F := F) l]
  ring

/-- ★★★★**`p₃ = (l−1)(3−l)/8`**。 -/
theorem muPow_three_zero {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    muPow l ζ 3 0 = ((l : F) - 1) * (3 - (l : F)) / 8 := by
  refine mul_left_cancel₀ (cast_prime_ne_zero (F := F) hl) ?_
  rw [muPow_zero_from_M hl hζ 2, muM_two_one hl hζ]
  ring

/-- `M_3^{(1)} = l(l−1)(l−2)(l−3)(9−l)/80`。 -/
theorem muM_one_three {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    muM l ζ 1 3
      = (l : F) * ((l : F) - 1) * ((l : F) - 2) * ((l : F) - 3) * (9 - (l : F)) / 80 := by
  rw [muM_rec hζ 0 3, muPow_one_zero hl hζ, muM_zero_four hl hζ,
    cast_choose_four_eq (F := F) l]
  ring

/-- `M_2^{(2)} = l(l−1)(l−2)(−l²−48l+193)/720`。 -/
theorem muM_two_two {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    muM l ζ 2 2
      = (l : F) * ((l : F) - 1) * ((l : F) - 2)
          * (-(l : F) ^ 2 - 48 * (l : F) + 193) / 720 := by
  rw [muM_rec hζ 1 2, muPow_two_zero hl hζ, muM_one_three hl hζ,
    cast_choose_three_eq (F := F) l]
  ring

/-- `M_1^{(3)} = l(l−1)(l³+l²−109l+251)/720`。 -/
theorem muM_three_one {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    muM l ζ 3 1
      = (l : F) * ((l : F) - 1)
          * ((l : F) ^ 3 + (l : F) ^ 2 - 109 * (l : F) + 251) / 720 := by
  rw [muM_rec hζ 2 1, muPow_three_zero hl hζ, muM_two_two hl hζ,
    cast_choose_two_eq (F := F) l]
  ring

/-- ★★★★★★★★★★**`p₄ = (l−1)(l³+l²−109l+251)/720`**——新しい種は要らなかった。 -/
theorem muPow_four_zero {l : ℕ} (hl : l.Prime) {ζ : F} (hζ : IsPrimitiveRoot ζ l) :
    muPow l ζ 4 0
      = ((l : F) - 1) * ((l : F) ^ 3 + (l : F) ^ 2 - 109 * (l : F) + 251) / 720 := by
  refine mul_left_cancel₀ (cast_prime_ne_zero (F := F) hl) ?_
  rw [muPow_zero_from_M hl hζ 3, muM_three_one hl hζ]
  ring

end CharZeroChain

/-! ## ★出典の紐付け(`.src`) -/

def muM.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(重み付きの和 M_m^{(k)}。★無条件)",
    sectionId := "genell-lemma-3-2" }

def muM_rec.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(層をひとつ上げる重み付きの漸化式。★ζ が原始 l 乗根)",
    sectionId := "genell-lemma-3-2" }

def muM_level_zero.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(M_m^{(0)} = l·C(l−1,m) − C(l,m+1)。★l は素数)",
    sectionId := "genell-lemma-3-2" }

def muPow_zero_from_M.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(l·A^{(k+1)}_0 = M_1^{(k)}。★l は素数)",
    sectionId := "genell-lemma-3-2" }

/-! ## ★出典の紐付け(`.src`) -/

def muPow.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(A^{(k)}_a = ∑_{ζ≠1} ζ^a/(1−ζ)^k。★無条件)",
    sectionId := "genell-lemma-3-2" }

def muPow_succ_sub.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(一段下げる。★ζ が原始 l 乗根)",
    sectionId := "genell-lemma-3-2" }

def muPow_sum_range.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(a について一周すると消える。★l は素数)",
    sectionId := "genell-lemma-3-2" }

def muPow_zero_rec.src : ABC3.Meta.Source :=
  { paper := "GenEll", pdfPage := 15,
    item := "Lemma 3.2, (ii)(層をひとつ上げる漸化式。★l は素数)",
    sectionId := "genell-lemma-3-2" }

end ABC3.Found.GaloisRep
